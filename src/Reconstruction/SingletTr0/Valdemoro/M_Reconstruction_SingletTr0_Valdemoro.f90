! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief SingletTr0 Valdemoro reconstruction scheme.
!>
!> @details
!> Implements the Valdemoro approximation for trace-free singlet states. Operates
!> on rdm2abTr0 and rdm1a separately, producing both collision operator components.
!>
!> **JSON path**: `reconstruction.singletTr0.valdemoro`
!>
!> The output consists of:
!> - collOp2ab: Trace-free collision contribution to d(rdm2abTr0)/dt
!> - collOp1a: Collision contribution to d(rdm1a)/dt
module M_Reconstruction_SingletTr0_Valdemoro
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the singlet Valdemoro implementation at runtime.
    !>
    !> Wires the singlet-specific procedure pointers (e.g., alpha–beta block
    !> builders) to Valdemoro implementations and may read/apply variant options
    !> as required by the codebase’s configuration system.
    module subroutine Reconstruction_SingletTr0_Valdemoro_Fabricate
    end subroutine
  end interface

end module
