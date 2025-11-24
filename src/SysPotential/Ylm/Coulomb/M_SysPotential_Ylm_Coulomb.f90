! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Central Coulomb potential in a Ylm representation.
!>
!> Parameterizes a 1/r potential of total charge `SysPotential_Ylm_Coulomb_charge`
!> and wires a Ylm-compatible implementation via the module factory.
module M_SysPotential_Ylm_Coulomb
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> set charge from config and wire the Ylm Coulomb implementation.
    module subroutine SysPotential_Ylm_Coulomb_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Total charge Z of the central Coulomb potential (V(r) = -Z / r in a.u.).
  real(R64) :: SysPotential_Ylm_Coulomb_charge

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
