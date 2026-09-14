! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Dispatcher for generic (spin-agnostic) purification backends.
!>
!> @details
!> When the JSON configuration contains `purification.generic`, this module's
!> `Purification_Generic_Fabricate` routine is invoked to select the specific
!> algorithm (currently only `alcoba` is available).
!>
!> @par Available Methods
!> | JSON Key                         | Backend                            | Status |
!> |----------------------------------|------------------------------------|--------|
!> | `purification.generic.alcoba`    | `M_Purification_Generic_Alcoba`    | Stub   |
!>
!> @par Use Case
!> The generic branch is intended for systems where spin symmetry is broken
!> or not exploited. The Alcoba method (when fully implemented) will operate
!> directly on the full unrestricted 2-RDM without spin-block decomposition.
module M_Purification_Generic
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Inspect `purification.generic.*` configuration and
    !> delegate to a concrete generic backend (e.g., Alcoba) which will assign
    !> the procedure pointers in `M_Purification`.
    module subroutine Purification_Generic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module

