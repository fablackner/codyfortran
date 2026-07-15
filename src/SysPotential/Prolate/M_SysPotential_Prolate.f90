! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Prolate-grid external potential backend (two-center coordinates).
!>
!> Provides the interface and state for potentials represented in azimuthal
!> channels up to `SysPotential_Prolate_mmax` on the prolate-spheroidal grid.
!> Concrete implementations (e.g., the two-center nuclear Coulomb attraction)
!> are selected and wired by `SysPotential_Prolate_Fabricate`.
module M_SysPotential_Prolate
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Prolate factory: parse config and wire concrete implementations.
    module subroutine SysPotential_Prolate_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Maximum azimuthal quantum number used to represent the potential
  !> (0 for phi-independent potentials such as the nuclear attraction).
  integer(I32) :: SysPotential_Prolate_mmax

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Compute the spatial component V_m(xi, eta; t) of the external potential
  !> in the e^(i m phi)/sqrt(2 pi) channel convention.
  procedure(I_SysPotential_Prolate_FillExternalPotentialChannel), pointer :: SysPotential_Prolate_FillExternalPotentialChannel
  abstract interface
    subroutine I_SysPotential_Prolate_FillExternalPotentialChannel(potM, m, time, bt_)
      import :: I32, R64
      !> Output: spatial channel array of size nXi*nEta (xi fastest).
      complex(R64), intent(out)          :: potM(:)
      !> Azimuthal quantum number (|m| <= SysPotential_Prolate_mmax).
      integer(I32), intent(in)           :: m
      !> Current simulation time.
      real(R64), intent(in)              :: time
      !> Optional body type selector.
      integer(I32), intent(in), optional :: bt_
    end subroutine
  end interface

end module
