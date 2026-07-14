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
!> Uses modified Gram–Schmidt with one re-orthogonalization pass ("twice is
!> enough") and metric-aware inner products:
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
!> Orthonormalizes the columns of fSet in place. Each column i is made
!> orthogonal to all preceding columns j < i by sequentially subtracting
!> projections (modified Gram–Schmidt), then normalized. A second projection
!> sweep re-orthogonalizes each vector before normalization: a single MGS
!> sweep leaves a residual overlap of order eps·kappa (condition number of
!> the input set), while the "twice is enough" strategy (Giraud et al.)
!> restores orthogonality to machine precision for any numerically
!> non-degenerate input. The metric is provided by the active back-end's
!> Grid_InnerProduct.
!>
!> @param[in,out] fSet  Matrix of size (Grid_nPoints, nSet) containing the
!>                      column vectors to orthonormalize.
  subroutine Orthonormalize(fSet)

    complex(R64), intent(inout), contiguous  :: fSet(:, :)

    integer(I32) :: i1, j1, iPass, nSet
    real(R64) :: norm
    complex(R64) :: ovlp

    nSet = size(fSet, 2)

    do i1 = 1, nSet

      ! Two modified Gram–Schmidt sweeps against all earlier vectors
      do iPass = 1, 2
        do j1 = 1, i1 - 1

          ovlp = Grid_InnerProduct(fSet(:, j1), fSet(:, i1))

          fSet(:, i1) = fSet(:, i1) - fSet(:, j1) * ovlp
        end do
      end do

      norm = real(Grid_InnerProduct(fSet(:, i1), fSet(:, i1)), kind=R64)

      if (.not. (norm > 0.0_R64)) then
        error stop "Grid Orthonormalize: vector with non-positive norm (linearly dependent input set?)"
      end if

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
    use M_Utils_BlasLib

    complex(R64), intent(inout), contiguous  :: fSetInOut(:, :)
    complex(R64), intent(in), contiguous     :: fSet(:, :)

    integer(I32) :: i1, j1, nSet, nG
    complex(R64), allocatable :: inp(:, :)

    nG = size(fSet, 1)
    nSet = size(fSet, 2)

    allocate (inp(nSet, nSet))

    !$omp parallel do default(shared) private(i1, j1) collapse(2)
    do j1 = 1, nSet
      do i1 = 1, nSet

        inp(i1, j1) = Grid_InnerProduct(fSet(:, i1), fSetInOut(:, j1))

      end do
    end do
    !$omp end parallel do

    ! Subtract projection onto each basis vector:
    !   fSetInOut = fSetInOut - fSet * inp   (single rank-nSet update via ZGEMM)
    ! WARNING: This procedure is delicate and tends to lose orthonormality
    ! over long propagations because numerical errors accumulate and interact
    ! with the nonlinear coupling between orbital and coefficient dynamics.
    call ZGEMM('N', 'N', nG, nSet, nSet, &
               cmplx(-1.0_R64, 0.0_R64, kind=R64), fSet, nG, inp, nSet, &
               cmplx(1.0_R64, 0.0_R64, kind=R64), fSetInOut, nG)

  end subroutine

end submodule
