! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Ylm.f90
!> @brief Spherical-harmonics kinetic operators (radial Laplacian + centrifugal).
!>
!> @details
!> This module covers kinetic operators in a **spherical-harmonics (Y_l^m)**
!> representation, where the spatial basis is factored as:
!>
!>     ψ(r,θ,φ) = Σ_{l,m} f_{lm}(r) · Y_l^m(θ,φ)
!>
!> The kinetic energy operator in spherical coordinates separates into:
!>
!>     T̂ = −(1/2m) [ (1/r²) ∂/∂r (r² ∂/∂r) − L̂²/(2m r²) ]
!>
!> where L̂² has eigenvalue l(l+1). For each (l,m) channel, this becomes:
!>
!>     T̂_{lm} f(r) = −(1/2m) [ d²f/dr² + (2/r) df/dr − l(l+1)/r² f ]
!>
!> Implementation Notes
!> --------------------
!> The radial part uses the standard transformation g(r) = r·f(r), which converts
!> the radial Laplacian to a simple second derivative:
!>
!>     [ d²/dr² + (2/r) d/dr ] f = (1/r) d²g/dr²
!>
!> This avoids the 2/r singularity and enables standard FD/FEDVR discretization.
!>
!> @see M_SysKinetic_Ylm_Laplacian, M_Grid_Ylm
module M_SysKinetic_Ylm
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Fabricate the Ylm kinetic backend and bind the radial operator.
    !>
    !> Currently supports: `laplacian` (with FinDiff or FEDVR sub-backends).
    module subroutine SysKinetic_Ylm_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to the radial kinetic operator for a single (l,m) channel.
  !>
  !> Applies: T̂_{lm} f(r) = −(1/2m) [ d²f/dr² + (2/r)df/dr − l(l+1)/r² f ]
  procedure(I_SysKinetic_Ylm_MultiplyWithRadialKineticOp), pointer :: SysKinetic_Ylm_MultiplyWithRadialKineticOp
  abstract interface
    !> @brief Apply the radial kinetic operator to a single (l,m) channel.
    !>
    !> @param[out] dOrbLm  Output radial function after applying T̂_{lm}
    !> @param[in]  orbLm   Input radial function f(r) for the (l,m) channel
    !> @param[in]  l       Angular momentum quantum number l ≥ 0
    !> @param[in]  m       Magnetic quantum number −l ≤ m ≤ l
    !> @param[in]  time    Current simulation time (reserved for extensions)
    !> @param[in]  bt_     Optional body-type index for mass scaling (default: 1)
    subroutine I_SysKinetic_Ylm_MultiplyWithRadialKineticOp(dOrbLm, orbLm, l, m, time, bt_)
      import :: I32, R64
      complex(R64), intent(out) :: dOrbLm(:)
      complex(R64), intent(in) :: orbLm(:)
      integer(I32), intent(in)  :: l
      integer(I32), intent(in)  :: m
      real(R64), intent(in) :: time
      integer(I32), intent(in), optional :: bt_
    end subroutine
  end interface

end module
