! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Spherical-ylm kinetic operators and setup.
!>
!> This module covers kinetic operators in a spherical-harmonics (Y_l^m)
!> representation, focusing on the radial part of the Laplacian. Procedure
!> pointers are assigned at runtime by `SysKinetic_Ylm_Fabricate`.
module M_SysKinetic_Ylm
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the ylm kinetic backend and bind the radial operator.
    module subroutine SysKinetic_Ylm_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the routine that applies the radial part of the Laplacian.
  procedure(I_SysKinetic_Ylm_MultiplyWithRadialKineticOp), pointer :: SysKinetic_Ylm_MultiplyWithRadialKineticOp
  abstract interface
    !> Applies the radial operator associated with −(1/2m)[∂²/∂r² + (2/r)∂/∂r − l(l+1)/r²]
    !> to the radial function of a single (l, m) channel.
    !>
    !> Implementation note: using the standard transformation g(r) = r f(r),
    !> the radial Laplacian acting on f can be represented via a second
    !> derivative on g, i.e. Laplacian(f) = (1/r) d²g/dr² − l(l+1) f / r².
    subroutine I_SysKinetic_Ylm_MultiplyWithRadialKineticOp(dOrbLm, orbLm, l, m, time, bt_)
      import :: I32, R64
      !> Output radial channel after applying the kinetic operator
      complex(R64), intent(out) :: dOrbLm(:)
      !> Input radial function f(r) for the (l, m) channel
      complex(R64), intent(in) :: orbLm(:)
      !> Angular momentum quantum number l
      integer(I32), intent(in)  :: l
      !> Magnetic quantum number m
      integer(I32), intent(in)  :: m
      !> Current time (supports time-dependent masses/masks if configured)
      real(R64), intent(in) :: time
      !> Optional body-type index for mass scaling (1/(2 m_bt))
      integer(I32), intent(in), optional :: bt_
    end subroutine
  end interface

end module
