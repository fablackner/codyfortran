! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ylm grid API for spherical-harmonic representations and runtime wiring.
!>
!> Provides a sparse/structured representation where functions are expressed on
!> a radial grid with spherical-harmonic (l,m) channels. The Fabricate routine
!> selects and initializes a concrete radial back-end (e.g., constant, FEDVR,
!> FEDVR-ECS) and wires operators acting in this mixed basis.
module M_Grid_Ylm
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a Ylm grid back-end.
    !>
    !> Initializes radial nodes/weights and angular layout up to lmax, and
    !> associates Ylm-specific procedures defined below.
    module subroutine Grid_Ylm_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Maximum angular momentum quantum number lmax (0 <= l <= lmax).
  integer(I32) :: Grid_Ylm_lmax

  !> Maximum radius value rmax of the radial grid.
  real(R64) :: Grid_Ylm_rmax

  !> Total integration weights for each (r,l,m) channel in flattened order.
  !> Combines radial weights and angular quadrature/normalization as used by
  !> the active Ylm back-end.
  real(R64), allocatable :: Grid_Ylm_weights(:)

  !> Component-wise coordinates for each (r,l,m) entry in flattened order.
  real(R64), allocatable :: Grid_Ylm_rCoord(:)
  integer(I32), allocatable :: Grid_Ylm_lCoord(:)
  integer(I32), allocatable :: Grid_Ylm_mCoord(:)

  !> Number of radial grid points.
  integer(I32) :: Grid_Ylm_nRadial

  !> Radial quadrature nodes r_i.
  real(R64), allocatable :: Grid_Ylm_radialPoints(:)

  !> Radial quadrature weights w_i^r.
  real(R64), allocatable :: Grid_Ylm_radialWeights(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the spatial-product routine producing a Ylm expansion.
  procedure(I_Grid_Ylm_SpatialProduct), pointer :: Grid_Ylm_SpatialProduct
  abstract interface
    !> Compute fOut = sum_r fConjg(r,theta,phi)^* f(r,theta,phi) expanded in Y_lm.
    !>
    !> Produces the spherical-harmonic coefficients up to lMaxOut from the
    !> product of two fields that may each be represented with their own lmax.
    !> Optional flags:
    !>  - conjgQ_: when true, treats fConjg as already conjugated (skips extra).
    !>  - withWeightsQ_: when true, includes the appropriate metric weights.
    pure subroutine I_Grid_Ylm_SpatialProduct(fOut, fConjg, f, lMaxOut, lMaxConjg, lMax, conjgQ_, withWeightsQ_)
      import :: I32, R64
      complex(R64), intent(out), contiguous :: fOut(:)
      complex(R64), intent(in), contiguous  :: fConjg(:)
      complex(R64), intent(in), contiguous  :: f(:)
      integer(I32), intent(in)              :: lMaxOut
      integer(I32), intent(in)              :: lMaxConjg
      integer(I32), intent(in)              :: lMax
      logical, intent(in), optional         :: conjgQ_
      logical, intent(in), optional         :: withWeightsQ_
    end subroutine
  end interface

  !> Pointer to the setup procedure for initializing the Ylm grid system.
  procedure(I_Grid_Ylm_Setup), pointer :: Grid_Ylm_Setup
  abstract interface
    !> Initialize the active Ylm back-end (radial grid and angular layout).
    subroutine I_Grid_Ylm_Setup
    end subroutine
  end interface

  !> Pointer to extract a specific (l,m) component from a function on the grid.
  procedure(I_Grid_Ylm_GetLmComponent), pointer :: Grid_Ylm_GetLmComponent
  abstract interface
    !> Extract the radial function f_lm(r) from a full Ylm-represented field.
    pure subroutine I_Grid_Ylm_GetLmComponent(fLmOut, l, m, f)
      import :: I32, R64
      !> Output radial function for the specified (l,m) component.
      complex(R64), intent(out) :: fLmOut(:)
      !> The angular momentum quantum number.
      integer(I32), intent(in) :: l
      !> The magnetic quantum number.
      integer(I32), intent(in) :: m
      !> The full function defined on the Ylm grid (flattened order).
      complex(R64), intent(in) :: f(:)
    end subroutine
  end interface

  !> Pointer to set/overwrite a specific (l,m) component in a function on the grid.
  procedure(I_Grid_Ylm_SetLmComponent), pointer :: Grid_Ylm_SetLmComponent
  abstract interface
    !> Set f_out so that its (l,m) channel equals the provided radial function.
    pure subroutine I_Grid_Ylm_SetLmComponent(fOut, l, m, fLm)
      import :: I32, R64
      !> The full function on the Ylm grid after setting the component.
      complex(R64), intent(out) :: fOut(:)
      !> The angular momentum quantum number.
      integer(I32), intent(in) :: l
      !> The magnetic quantum number.
      integer(I32), intent(in) :: m
      !> The input radial function for the specified (l,m) component.
      complex(R64), intent(in) :: fLm(:)
    end subroutine
  end interface

  !> Pointer to add/accumulate a specific (l,m) component into a function on the grid.
  procedure(I_Grid_Ylm_AddLmComponent), pointer :: Grid_Ylm_AddLmComponent
  abstract interface
    !> Add the provided radial function into the (l,m) channel of fInOut.
    pure subroutine I_Grid_Ylm_AddLmComponent(fInOut, l, m, fLm)
      import :: I32, R64
      !> The full function on the Ylm grid to which the component is added.
      complex(R64), intent(inout) :: fInOut(:)
      !> The angular momentum quantum number.
      integer(I32), intent(in) :: l
      !> The magnetic quantum number.
      integer(I32), intent(in) :: m
      !> The input radial function for the specified (l,m) component.
      complex(R64), intent(in) :: fLm(:)
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
