! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Single-particle quantum dynamics (no mean-field or many-body terms).
!>
!> Evolves a single wavefunction under an external potential; no exchange or
!> correlation effects are included. Useful for testing grid implementations and
!> studying fundamental phenomena (wave packets, tunneling, scattering).
!>
!> State layout
!> - The packed state stores the complex single-particle wavefunction in the
!>   chosen representation (grid or basis), with layout determined during setup.
module M_Method_Sb
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds Sb procedure pointers at runtime and initializes the
    !> single-particle state. Typically sets `Method_Setup`, `Method_GetEnergy`,
    !> and `Method_TimeDerivative`.
    module subroutine Method_Sb_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
