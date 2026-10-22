! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Submodule S_CoeffsInit implements the fabrication logic for CI coefficient
!> initialization.
!>
!> Branching Logic
!> ---------------
!> The factory dispatches to exactly one concrete initializer based on
!> JSON configuration presence (checked in priority order):
!>
!>   1. `coeffsInit.load`    → Load coefficients from binary file
!>   2. `coeffsInit.unary`   → Single-configuration ground state (c₁ = 1)
!>   3. `coeffsInit.excited` → Excited state via creation/annihilation ops
!>
!> Mutual exclusivity is enforced: only the first matching key is honored.
!> If no key is found, fabrication fails with an error stop.
!>
!> `CoeffsInit_Initialize` is bound to the local `Initialize`, which calls the
!> concrete initializer and then projects on the declared total-M sector
!> (Coeffs_ProjectOnSector, NoOp without coeffs.generic.totalM). This is the
!> CI counterpart of Orbs_ProjectOnSymmetry in OrbsInit_Initialize.
submodule(M_CoeffsInit) S_CoeffsInit

  implicit none

  !> Concrete initializer bound by the delegated fabricator
  procedure(I_CoeffsInit_Initialize), pointer, save :: InitializeConcrete => null()

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Inspects JSON configuration and delegates to the appropriate concrete
  !> initializer's fabrication routine.
  !>
  !> The procedure pointer `CoeffsInit_Initialize` (and optionally
  !> `CoeffsInit_Setup`) in M_CoeffsInit will be bound by the delegated
  !> fabricator.
  !>
  !> @note The branching order matters: `load` takes precedence over `unary`,
  !>       which takes precedence over `excited`. This allows fallback
  !>       configurations where multiple keys might coexist.
  module subroutine CoeffsInit_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_CoeffsInit_Load
    use M_CoeffsInit_Unary
    use M_CoeffsInit_Excited

    call Say_Fabricate("coeffsInit")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch: dispatch to concrete initializer
    !------------------------------------

    if (Json_GetExistence("coeffsInit.load")) then
      call CoeffsInit_Load_Fabricate

    else if (Json_GetExistence("coeffsInit.unary")) then
      call CoeffsInit_Unary_Fabricate

    else if (Json_GetExistence("coeffsInit.excited")) then
      call CoeffsInit_Excited_Fabricate

    else
      error stop "coeffsInit is missing one of: load, unary, excited"
    end if

    InitializeConcrete => CoeffsInit_Initialize
    CoeffsInit_Initialize => Initialize

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Concrete initializer followed by the projection on the declared total-M
  !> sector; stops if the initial state does not lie in it (e.g. a unary or
  !> excited configuration of other total M, or a coeffs.in of another sector).
  subroutine Initialize(coeffs)
    use M_Coeffs

    complex(R64), intent(out), contiguous :: coeffs(:)

    call InitializeConcrete(coeffs)
    call Coeffs_ProjectOnSector(coeffs)

  end subroutine

end submodule
