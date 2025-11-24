! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default TDHx wiring for the ground-state solver.
!>
!> This module provides the factory routine that binds TDHx callbacks to a
!> reference/default implementation (e.g., CPU, baseline linear algebra).
!> It contains no heavy numerical kernels itself; it only assigns the
!> procedure pointers exported by the TDHx interface modules.
module M_GroundSolver_Tdhx_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind TDHx procedure pointers to the default implementation.
    !>
    !> Responsibilities typically include assigning the Hartree–Fock action
    !> routine and any auxiliary callbacks needed by the ground-state solver.
    module subroutine GroundSolver_Tdhx_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module

