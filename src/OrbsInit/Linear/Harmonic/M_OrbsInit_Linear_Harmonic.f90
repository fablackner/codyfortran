! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module M_OrbsInit_Linear_Harmonic
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate a harmonic-oscillator based linear initializer.
    !> Configures parameters for harmonic orbitals on a 1D grid and wires the
    !> callouts in the linear backend (e.g., ground/excited states centered at a
    !> specified position with frequency `omega`).
    module subroutine OrbsInit_Linear_Harmonic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Center position x0 of the harmonic potential (same units as grid coordinate).
  real(R64) :: OrbsInit_Linear_Harmonic_position
  !> Angular frequency ω of the harmonic potential (units consistent with model).
  real(R64) :: OrbsInit_Linear_Harmonic_omega

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
