! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Factory for split-operator propagator backends.
!>
!> Dispatches to the appropriate composition order:
!> - `order2`: Symmetric Strang splitting, O(Δt³) local error
!> - `order4`: Yoshida/Forest–Ruth composition, O(Δt⁵) local error
!>
!> The split-step method decomposes Ĥ = Â + B̂ and approximates:
!>   exp(−iĤΔt) ≈ product of exp(−iÂ·aₖΔt) and exp(−iB̂·bₖΔt)
!>
!> Each operator (Â, B̂) is applied via a separate integrator in IntegratorList.
submodule(M_Propagator_SplitStep) S_Propagator_SplitStep

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Reads `"propagator.splitStep"` and delegates to order-specific backend.
!>
!> Exactly one order key must be present:
!>   - `"propagator.splitStep.order2"` → 2nd-order Strang splitting
!>   - `"propagator.splitStep.order4"` → 4th-order Yoshida composition
  module subroutine Propagator_SplitStep_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Propagator_SplitStep_Order2
    use M_Propagator_SplitStep_Order4

    call Say_Fabricate("propagator.splitStep")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch: dispatch to order-specific backend
    !------------------------------------

    if (Json_GetExistence("propagator.splitStep.order2")) then
      call Propagator_SplitStep_Order2_Fabricate

    else if (Json_GetExistence("propagator.splitStep.order4")) then
      call Propagator_SplitStep_Order4_Fabricate

    else
      error stop "propagator.splitStep is missing one of: order2, order4"
    end if

  end subroutine

end submodule
