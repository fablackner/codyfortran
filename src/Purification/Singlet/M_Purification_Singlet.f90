! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Dispatcher for singlet spin-symmetry purification backends.
!>
!> @details
!> When the JSON configuration contains `purification.singlet`, this module's
!> `Purification_Singlet_Fabricate` routine is invoked to select the specific
!> singlet algorithm (currently only `reduce` is available).
!>
!> @par Available Methods
!> | JSON Key                         | Backend                            |
!> |----------------------------------|------------------------------------|
!> | `purification.singlet.reduce`    | `M_Purification_Singlet_Reduce`    |
!>
!> @par Singlet Symmetry Context
!> Singlet (S=0) spin eigenstates have a 2-RDM that decomposes into symmetric
!> (singlet) and antisymmetric (triplet) blocks under particle exchange. This
!> structure is exploited by purification algorithms to separately enforce
!> positivity in each spin sector.
module M_Purification_Singlet
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Inspect the singlet branch of the configuration and delegate to a
    !> concrete singlet backend (e.g., Reduce) which will assign the procedure
    !> pointers in `M_Purification`.
    module subroutine Purification_Singlet_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module

