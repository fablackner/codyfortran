! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Generic grid API and runtime wiring for spatial representations.
!>
!> M_Grid provides a thin, implementation-agnostic interface for spatial
!> grids used throughout the code base. It exposes a small set of procedure
!> pointers (setup, inner product, orthonormalization, projection) and a
!> single size scalar. Concrete grid back-ends (1D constant, 2D/3D regular,
!> spherical/Ylm, FEDVR, …) are selected at runtime by Grid_Fabricate which
!> associates these procedure pointers and initializes module state owned by
!> the selected back-end.
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
