! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the root purification fabrication logic.
!>
!> @details
!> This submodule inspects the JSON configuration to determine which
!> purification family to activate and delegates to the corresponding
!> branch fabricator. The dispatch hierarchy is:
!>
!> @verbatim
!>   Purification_Fabricate
!>       │
!>       ├── purification.singlet    → Purification_Singlet_Fabricate
!>       │                                   └── reduce → Singlet_Reduce
!>       │
!>       ├── purification.singletTr0 → Purification_SingletTr0_Fabricate
!>       │                                   └── reduce → SingletTr0_Reduce
!>       │
!>       └── purification.generic    → Purification_Generic_Fabricate
!>                                           └── alcoba → Generic_Alcoba
!> @endverbatim
!>
!> If no `purification.*` key is found, a diagnostic is printed and
!> purification remains disabled (no-op).
submodule(M_Purification) S_Purification

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Dispatch to the appropriate purification family based on JSON config.
  module subroutine Purification_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Purification_Singlet
    use M_Purification_Generic
    use M_Purification_SingletTr0

    call Say_Fabricate("purification")

    !------------------------------------
    ! branch: select purification family
    !------------------------------------

    if (Json_GetExistence("purification.singlet")) then
      ! Singlet spin-symmetry adapted purification
      call Purification_Singlet_Fabricate

    else if (Json_GetExistence("purification.generic")) then
      ! Spin-agnostic (generic) purification
      call Purification_Generic_Fabricate

    else if (Json_GetExistence("purification.singletTr0")) then
      ! SingletTr0: singlet with explicit trace-free 2-RDM representation
      call Purification_SingletTr0_Fabricate

    else
      ! No purification configured - pointers remain bound to no-ops
      write (*, '(1X, A, A, A)') red, 'de: Purification => noPurification', reset

    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule

