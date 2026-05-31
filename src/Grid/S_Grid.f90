! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_Grid.f90
!> @brief Implementation submodule for M_Grid: fabrication and shared routines.
!>
!> @details
!> Implements Grid_Fabricate (the top-level factory that dispatches to the
!> appropriate back-end) and provides generic Gram–Schmidt orthonormalization
!> and subspace projection routines that are shared across all coordinate
!> back-ends.
!>
!> @par Fabrication Dispatch
!> The JSON configuration is probed in the following order:
!>   1. `grid.linear`    → Grid_Linear_Fabricate
!>   2. `grid.square`    → Grid_Square_Fabricate
!>   3. `grid.polar`     → Grid_Polar_Fabricate
!>   4. `grid.spherical` → Grid_Spherical_Fabricate
!>   5. `grid.ylm`       → Grid_Ylm_Fabricate
!>   6. `grid.lattice`   → Grid_Lattice_Fabricate
!>
!> @par Orthonormalization Algorithm
!> Uses the classical Gram–Schmidt procedure with metric-aware inner products:
!> @f[
!>   |\phi_i\rangle \leftarrow |\phi_i\rangle
!>     - \sum_{j<i} \langle\phi_j|\phi_i\rangle |\phi_j\rangle, \quad
!>   |\phi_i\rangle \leftarrow \frac{|\phi_i\rangle}{\sqrt{\langle\phi_i|\phi_i\rangle}}
!> @f]
!>
!> @par Subspace Projection
!> Removes the component of each input vector along a given orthonormal basis:
!> @f[
!>   |\psi_j\rangle \leftarrow |\psi_j\rangle
!>     - \sum_{i} \langle\phi_i|\psi_j\rangle |\phi_i\rangle
!> @f]
!> This operation is delicate in time propagation because numerical errors
!> can accumulate and violate gauge constraints; see comments in the code.
submodule(M_Grid) S_Grid

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Top-level factory: parse JSON and dispatch to the appropriate back-end.
!>
!> Binds the generic Grid_Orthonormalize and Grid_ProjectOnSubspace pointers
!> to the shared implementations in this submodule, then delegates to the
!> selected coordinate family (Linear, Square, Polar, Spherical, Ylm, Lattice).
  module subroutine Grid_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid_Linear
    use M_Grid_Square
    use M_Grid_Polar
    use M_Grid_Spherical
    use M_Grid_Ylm
    use M_Grid_Lattice

    call Say_Fabricate("grid")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Grid_Orthonormalize => Orthonormalize
    Grid_ProjectOnSubspace => ProjectOnSubspace

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("grid.linear")) then
      call Grid_Linear_Fabricate

    else if (Json_GetExistence("grid.square")) then
      call Grid_Square_Fabricate

    else if (Json_GetExistence("grid.polar")) then
      call Grid_Polar_Fabricate

    else if (Json_GetExistence("grid.spherical")) then
      call Grid_Spherical_Fabricate

    else if (Json_GetExistence("grid.ylm")) then
      call Grid_Ylm_Fabricate

    else if (Json_GetExistence("grid.lattice")) then
      call Grid_Lattice_Fabricate

    else
      error stop "grid is missing one of: linear, square, polar, spherical, ylm, lattice"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Modified Gram–Schmidt orthonormalization using Grid_InnerProduct.
!>
!> @details
!> Orthonormalizes the columns of fSet in place. Each column i is first made
!> orthogonal to all preceding columns j < i by subtracting projections, then
!> normalized to unit length. The metric is provided by the active back-end's
!> Grid_InnerProduct.
!>
!> @warning For large sets or poorly conditioned bases, classical Gram–Schmidt
!> can lose orthogonality. Consider re-orthogonalization or QR-based methods
!> if numerical stability is critical.
!>
!> @param[in,out] fSet  Matrix of size (Grid_nPoints, nSet) containing the
!>                      column vectors to orthonormalize.
  subroutine Orthonormalize(fSet)

    complex(R64), intent(inout), contiguous  :: fSet(:, :)

    integer(I32) :: i1, j1, nSet
    real(R64) :: norm
    complex(R64) :: ovlp

    nSet = size(fSet, 2)

    ! Classical Gram–Schmidt orthonormalization
    ! For each vector i: subtract projections onto all earlier vectors j,
    ! then normalize to unit length.
    do i1 = 1, nSet
      do j1 = 1, i1 - 1

        ovlp = Grid_InnerProduct(fSet(:, j1), fSet(:, i1))

        fSet(:, i1) = fSet(:, i1) - fSet(:, j1) * ovlp
      end do

      norm = real(Grid_InnerProduct(fSet(:, i1), fSet(:, i1)), kind=R64)
      fSet(:, i1) = fSet(:, i1) / sqrt(norm)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Project vectors onto the orthogonal complement of a given subspace.
!>
!> @details
!> For each column in fSetInOut, removes its component along the subspace
!> spanned by the columns of fSet. Assumes fSet is already orthonormal w.r.t.
!> Grid_InnerProduct; if not, results are undefined.
!>
!> @note This operation is particularly delicate in MCTDH-style propagation
!> where gauge constraints require the time-derivative orbitals to remain
!> orthogonal to the occupied subspace. Accumulated numerical errors can break
!> this balance, especially when the correlation operator introduces
!> nonlinearities coupling orbital and coefficient dynamics.
!>
!> @param[in,out] fSetInOut  Vectors to project (Grid_nPoints x N).
!> @param[in]     fSet       Orthonormal basis spanning the subspace to
!>                           project out (Grid_nPoints x M).
  subroutine ProjectOnSubspace(fSetInOut, fSet)

    complex(R64), intent(inout), contiguous  :: fSetInOut(:, :)
    complex(R64), intent(in), contiguous     :: fSet(:, :)

    integer(I32) :: i1, j1, nSet
    complex(R64), allocatable :: inp(:, :)

    nSet = size(fSet, 2)

    allocate (inp(nSet, nSet))

    do j1 = 1, nSet
      do i1 = 1, nSet

        inp(i1, j1) = Grid_InnerProduct(fSet(:, i1), fSetInOut(:, j1))

      end do
    end do

    ! Subtract projection onto each basis vector.
    ! WARNING: This procedure is delicate and tends to lose orthonormality
    ! over long propagations because numerical errors accumulate and interact
    ! with the nonlinear coupling between orbital and coefficient dynamics.
    do j1 = 1, nSet
      do i1 = 1, nSet

        fSetInOut(:, j1) = fSetInOut(:, j1) - inp(i1, j1) * fSet(:, i1)

      end do
    end do

  end subroutine

end submodule
