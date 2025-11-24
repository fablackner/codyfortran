! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Linear-grid external potential backend.
!>
!> Exposes the linear-grid specific setup hook and a factory that wires concrete
!> external potential implementations (e.g., harmonic, soft-Yukawa) at runtime.
!> The selected model and parameters are taken from the JSON configuration.
module M_SysPotential_Linear
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Linear factory: parse config and wire concrete implementations.
    module subroutine SysPotential_Linear_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Linear-grid setup hook (allocations, masks, caches for linear grids).
  procedure(I_SysPotential_Linear_Setup), pointer :: SysPotential_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize the linear potential backend.
    subroutine I_SysPotential_Linear_Setup
    end subroutine
  end interface

end module
