! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Factory implementation for the Propagator facade.
!>
!> Reads the JSON configuration and dispatches to the appropriate backend:
!> - `single`: Direct time integration via IntegratorList (RK, SIL, etc.)
!> - `splitStep`: Suzuki–Trotter factorization (order2, order4)
!> - `eigenExpansion`: Diagonalization + phase evolution
!>
!> Exactly one backend key must be present under `"propagator"`.
submodule(M_Propagator) S_Propagator

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Reads the JSON `"propagator"` block and delegates to the matching backend.
!>
!> Backend selection is mutually exclusive. The first matching key wins:
!>   1. `"propagator.single"` → Propagator_Single_Fabricate
!>   2. `"propagator.splitStep"` → Propagator_SplitStep_Fabricate
!>   3. `"propagator.eigenExpansion"` → Propagator_EigenExpansion_Fabricate
!>
!> On success, the module-level procedure pointers `Propagator_Setup` and
!> `Propagator_Propagate` are bound to the chosen backend's implementations.
  module subroutine Propagator_Fabricate()
    use M_Utils_Json
    use M_Utils_Say
    use M_Propagator_Single
    use M_Propagator_SplitStep
    use M_Propagator_EigenExpansion

    call Say_Fabricate("propagator")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch: dispatch to backend based on JSON key
    !------------------------------------

    if (Json_GetExistence("propagator.single")) then
      call Propagator_Single_Fabricate

    else if (Json_GetExistence("propagator.splitStep")) then
      call Propagator_SplitStep_Fabricate

    else if (Json_GetExistence("propagator.eigenExpansion")) then
      call Propagator_EigenExpansion_Fabricate

    else
      error stop "propagator is missing one of: single, splitStep, eigenExpansion"
    end if

  end subroutine

end submodule
