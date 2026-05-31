! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_Grid.f90
!> @brief Generic grid API and runtime wiring for spatial representations.
!>
!> @details
!> M_Grid provides a thin, implementation-agnostic interface for spatial
!> discretizations used throughout the quantum many-body simulation framework.
!> It exposes:
!>   - A global scalar `Grid_nPoints` (total degrees of freedom)
!>   - Procedure pointers for setup, inner product, Gram–Schmidt
!>     orthonormalization, and subspace projection
!>
!> Concrete spatial back-ends are selected at runtime via `Grid_Fabricate`,
!> which reads the JSON configuration and delegates to the appropriate
!> implementation (Linear, Square, Polar, Spherical, Ylm, Lattice). Each
!> back-end defines its own metric/quadrature so that inner products,
!> orthonormality, and projection are metric-aware.
!>
!> @par Supported Back-ends
!>   | JSON Key         | Module             | Geometry / Basis                       |
!>   |------------------|--------------------|----------------------------------------|
!>   | `grid.linear`    | M_Grid_Linear      | 1D Cartesian (Const, FEDVR)            |
!>   | `grid.square`    | M_Grid_Square      | 2D Cartesian (Const)                   |
!>   | `grid.polar`     | M_Grid_Polar       | 2D polar (Const)                       |
!>   | `grid.spherical` | M_Grid_Spherical   | 3D spherical (Const)                   |
!>   | `grid.ylm`       | M_Grid_Ylm         | Radial + spherical harmonics Y_lm      |
!>   | `grid.lattice`   | M_Grid_Lattice     | 3D discrete lattice (Hubbard-style)    |
!>
!> @par Typical Initialization
!> @code{.f90}
!>   call Grid_Fabricate   ! selects & wires back-end from JSON
!>   call Grid_Setup       ! allocates coordinates, weights, etc.
!> @endcode
!>
!> @note Only import M_Grid in client code; implementation details live in
!>       the corresponding M_Grid_* and S_Grid_* modules.
!>
!> @see S_Grid.f90 for fabrication logic and default orthonormalization.
module M_Grid
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and wire a concrete grid back-end.
    !>
    !> Associates the procedure pointers below with the selected grid
    !> representation and initializes any required global state for later use
    !> (sizes, quadrature/metric, coordinate arrays, …). The active back-end
    !> defines the metric used by Grid_InnerProduct.
    module subroutine Grid_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Total number of degrees of freedom in the active spatial representation.
  !> For point-based grids this equals the number of grid points; for composite
  !> representations it is the flattened size used by the high-level solvers.
  integer(I32) :: Grid_nPoints

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup procedure for initializing the grid system.
  !> Must be associated by Grid_Fabricate before use.
  procedure(I_Grid_Setup), pointer :: Grid_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize the active grid back-end.
    !>
    !> Prepares module-level state owned by the back-end (coordinates,
    !> weights/metric, sizes, workspace). Safe to call more than once; a
    !> back-end may choose to reinitialize if configuration changed.
    subroutine I_Grid_Setup
    end subroutine
  end interface

  !> Pointer to the procedure for calculating the inner product of fields.
  !> The precise metric (e.g., quadrature weights, volume factors) is defined
  !> by the active back-end selected via Grid_Fabricate.
  procedure(I_Grid_InnerProduct), pointer :: Grid_InnerProduct
  abstract interface
    !> Calculate the inner product of two fields on the active grid.
    !>
    !> The result reflects the back-end’s metric. Arguments must be 1D views
    !> matching Grid_nPoints and be contiguous in memory.
    pure function I_Grid_InnerProduct(fConjg, f) result(res)
      import :: R64
      !> Inner product value according to the active grid metric.
      complex(R64) :: res
      !> Complex-conjugated first field (Dirac bra component).
      complex(R64), intent(in), contiguous :: fConjg(:)
      !> Second field (Dirac ket component).
      complex(R64), intent(in), contiguous :: f(:)
    end function
  end interface

  !> Pointer to orthonormalization routine.
  !> Uses the active grid metric as defined by the back-end.
  procedure(I_Grid_Orthonormalize), pointer :: Grid_Orthonormalize
  abstract interface
    !> Orthonormalize a set of column vectors in-place.
    !>
    !> fSet contains N columns of length Grid_nPoints. The routine enforces
    !> orthonormality with respect to Grid_InnerProduct and may use modified
    !> Gram–Schmidt or a back-end optimized variant.
    subroutine I_Grid_Orthonormalize(fSet)
      import :: R64
      !> Set of fields to be orthonormalized (size: Grid_nPoints x N).
      complex(R64), intent(inout), contiguous :: fSet(:, :)
    end subroutine
  end interface

  !> Pointer to projection-on-subspace routine.
  !> Projects using the active grid metric.
  procedure(I_Grid_ProjectOnSubspace), pointer :: Grid_ProjectOnSubspace
  abstract interface
    !> Project fSetInOut onto the orthogonal complement of span(fSet).
    !>
    !> Removes components of each column in fSetInOut along the subspace
    !> spanned by the columns of fSet, with orthogonality defined by
    !> Grid_InnerProduct.
    subroutine I_Grid_ProjectOnSubspace(fSetInOut, fSet)
      import :: R64
      !> Vectors to be projected (size: Grid_nPoints x N).
      complex(R64), intent(inout), contiguous :: fSetInOut(:, :)
      !> Basis defining the subspace (size: Grid_nPoints x M).
      complex(R64), intent(in), contiguous :: fSet(:, :)
    end subroutine
  end interface

end module
