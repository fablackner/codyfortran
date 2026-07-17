! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Spherical-harmonics gauge coupling: radial building blocks.
!>
!> In a Ylm representation ψ(r,θ,φ) = Σ_{lm} f_{lm}(r) Y_l^m(θ,φ), a
!> z-polarized vector potential couples neighboring angular momenta through
!>
!>     ∂_z [f_l Y_{lm}] = c_{l,m} (∂_r − l/r) f_l · Y_{l+1,m}
!>                      + c_{l−1,m} (∂_r + (l+1)/r) f_l · Y_{l−1,m}
!>
!> with c_{l,m} = sqrt( ((l+1)² − m²) / ((2l+1)(2l+3)) ). The angular
!> coupling structure lives in the gauge model (e.g., Velocity); this module
!> exports the radial building blocks bound by the discretization backends
!> (FEDVR, FEDVR-ECS): the radial first derivative and the (possibly complex)
!> radial coordinates entering the ±k/r terms.
module M_SysGauge_Ylm
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the Ylm gauge backend and dispatch to a gauge model.
    !>
    !> Currently supports: `velocity` (with Fedvr or FedvrEcs sub-backends).
    module subroutine SysGauge_Ylm_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Radial coordinates entering the ±k/r channel-coupling terms.
  !> Real radial points for Hermitian backends, complex contour points z(r)
  !> for the ECS contour. Filled by the active backend during SysGauge_Setup.
  complex(R64), allocatable :: SysGauge_Ylm_radialCoordinates(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Apply the radial first derivative to a single (l,m) channel:
  !> dfLm = d fLm / dr (on the ECS contour: d/dz along the contour).
  procedure(I_SysGauge_Ylm_ApplyRadialFirstDerivative), pointer :: SysGauge_Ylm_ApplyRadialFirstDerivative
  abstract interface
    !> Compute the radial first derivative of a radial channel function.
    subroutine I_SysGauge_Ylm_ApplyRadialFirstDerivative(dfLm, fLm)
      import :: R64
      !> Output: radial derivative of the channel function.
      complex(R64), intent(out) :: dfLm(:)
      !> Input: radial channel function f_{lm}(r).
      complex(R64), intent(in) :: fLm(:)
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
