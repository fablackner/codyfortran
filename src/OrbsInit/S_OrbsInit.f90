! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for M_OrbsInit.
!>
!> @details
!> Contains the fabrication logic that inspects JSON configuration and delegates
!> to the appropriate backend module (Linear, Lattice, GridPoint, Load, or Ylm).
!> Also provides the generic `Initialize` routine that loops over body types
!> and orbitals, calling the backend-specific `InitializeOrb` for each.
submodule(M_OrbsInit) S_OrbsInit

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Select and wire the orbital initialization backend from JSON config.
!>
!> @details
!> Inspects the JSON tree under "orbsInit" to determine which backend to use.
!> The search order is: linear → lattice → gridPoint → load → ylm.
!> Each backend's Fabricate routine further binds `OrbsInit_InitializeOrb`.
!>
!> JSON paths checked:
!> - orbsInit.linear    → OrbsInit_Linear_Fabricate
!> - orbsInit.lattice   → OrbsInit_Lattice_Fabricate
!> - orbsInit.gridPoint → OrbsInit_GridPoint_Fabricate
!> - orbsInit.load      → OrbsInit_Load_Fabricate
!> - orbsInit.ylm       → OrbsInit_Ylm_Fabricate
!> - orbsInit.prolate   → OrbsInit_Prolate_Fabricate
  module subroutine OrbsInit_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit_Linear
    use M_OrbsInit_GridPoint
    use M_OrbsInit_Lattice
    use M_OrbsInit_Load
    use M_OrbsInit_Ylm
    use M_OrbsInit_Prolate

    call Say_Fabricate("orbsInit")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_Initialize => Initialize

    !------------------------------------
    ! branch based on JSON configuration
    !------------------------------------

    if (Json_GetExistence("orbsInit.linear")) then
      call OrbsInit_Linear_Fabricate

    else if (Json_GetExistence("orbsInit.lattice")) then
      call OrbsInit_Lattice_Fabricate

    else if (Json_GetExistence("orbsInit.gridPoint")) then
      call OrbsInit_GridPoint_Fabricate

    else if (Json_GetExistence("orbsInit.load")) then
      call OrbsInit_Load_Fabricate

    else if (Json_GetExistence("orbsInit.ylm")) then
      call OrbsInit_Ylm_Fabricate

    else if (Json_GetExistence("orbsInit.prolate")) then
      call OrbsInit_Prolate_Fabricate

    else
      error stop "orbsInit is missing one of: linear, lattice, gridPoint, load, ylm, prolate"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Initialize all orbitals by looping over body types and orbital indices.
!>
!> @details
!> Iterates over all body types (e.g., spin species) and their respective orbital
!> counts. For each (body-type, orbital-index) pair, calls the backend-specific
!> `OrbsInit_InitializeOrb` to fill the corresponding column in the orbs array.
!>
!> The orbitals are stored column-major with layout:
!>   orbs(:, 1:nOrbs(1))                 → body type 1
!>   orbs(:, nOrbs(1)+1:nOrbs(1)+nOrbs(2)) → body type 2
!>   etc.
!>
!> @param[out] orbs  2D array (nGrid × nTotalOrbs) to be filled with orbital data.
!>                   Each column is a normalized orbital sampled on the grid.
  subroutine Initialize(orbs)
    use M_Utils_DataStorage
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Grid
    use M_Orbs

    complex(R64), intent(out), contiguous  :: orbs(:, :)

    integer(I32) :: ind  ! orbital index within current body type
    integer(I32) :: ibt  ! body type index
    integer(I32) :: i1   ! global orbital column index
    integer(I32) :: nBt  ! number of body types to initialize

    ! Restricted: only the shared spatial set (body type 1) is stored
    nBt = Method_Mb_nBodyTypes
    if (Orbs_restrictedQ) nBt = 1

    i1 = 0
    do ibt = 1, nBt
      do ind = 1, Method_Mb_OrbBased_nOrbs(ibt)
        i1 = i1 + 1

        call OrbsInit_InitializeOrb(orbs(:, i1), ind, ibt)

      end do
    end do

    ! NOTE: Orthonormalization is typically handled by Orbs module after init
    !call Orbs_Orthonormalize(orbs)

  end subroutine

end submodule
