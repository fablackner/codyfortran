! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Dispatcher for SingletTr0 spin-symmetry purification backends.
!>
!> @details
!> When the JSON configuration contains `purification.singletTr0`, this module's
!> `Purification_SingletTr0_Fabricate` routine is invoked to select the specific
!> algorithm (currently only `reduce` is available).
!>
!> @par Available Methods
!> | JSON Key                            | Backend                               |
!> |-------------------------------------|---------------------------------------|
!> | `purification.singletTr0.reduce`    | `M_Purification_SingletTr0_Reduce`    |
!>
!> @par SingletTr0 Representation
!> The SingletTr0 state container stores only the trace-free (cumulant) part
!> of the 2-RDM alongside an explicit 1-RDM. This reduces degrees of freedom
!> and is well-suited for systems where the 1-RDM is independently tracked.
module M_Purification_SingletTr0
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Inspect the singletTr0 branch of the configuration and delegate to a
    !> concrete singletTr0 backend (e.g., Reduce) which will assign the procedure
    !> pointers in `M_Purification`.
    module subroutine Purification_SingletTr0_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module

