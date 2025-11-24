! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> On-site interaction model for lattice grids.
!>
!> Provides the parameterization and fabrication hook for a local (delta-like)
!> on-site interaction used in lattice-based simulations.
module M_SysInteraction_Lattice_OnSite
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Register the on-site lattice interaction with the core system. Reads and
    !> validates parameters, then wires the concrete implementation.
    module subroutine SysInteraction_Lattice_OnSite_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Coupling constant of the local on-site interaction. Positive values
  !> indicate repulsion, negative values attraction (convention depends on
  !> the chosen Hamiltonian).
  real(R64) :: SysInteraction_Lattice_OnSite_Strength

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
