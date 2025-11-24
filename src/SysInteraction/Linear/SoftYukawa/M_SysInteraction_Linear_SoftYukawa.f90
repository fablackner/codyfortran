! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Softened, damped Yukawa interaction on a linear grid.
!>
!> Parameterization of a screened Coulomb/Yukawa-like kernel with optional
!> softening near the origin and exponential damping at large distances.
module M_SysInteraction_Linear_SoftYukawa
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Register the SoftYukawa model and expose its kernel to the linear-grid
    !> back-end (either direct or FFT-based).
    module subroutine SysInteraction_Linear_SoftYukawa_Fabricate
    end subroutine

    !> Kernel function K(r) of the SoftYukawa model evaluated at a scalar
    !> distance `r` (in the active grid's length units). The exact form is
    !> controlled by the parameters below.
    pure module function SysInteraction_Linear_SoftYukawa_Interaction(distance) result(res)
      real(R64) :: res
      real(R64), intent(in) :: distance
    end function
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Overall coupling strength of the kernel.
  real(R64) :: SysInteraction_Linear_SoftYukawa_strength
  !> Inner softening parameter (near-field regularization).
  real(R64) :: SysInteraction_Linear_SoftYukawa_softening1
  !> Outer softening/shape parameter.
  real(R64) :: SysInteraction_Linear_SoftYukawa_softening2
  !> Exponential damping coefficient controlling long-range decay.
  real(R64) :: SysInteraction_Linear_SoftYukawa_dampening

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
