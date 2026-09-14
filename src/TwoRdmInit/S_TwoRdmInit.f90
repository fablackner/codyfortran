! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_TwoRdmInit.f90
!> @brief Implementation submodule for M_TwoRdmInit.
!>
!> @details
!> This submodule provides:
!>   - `TwoRdmInit_Fabricate`: Factory dispatcher that reads JSON and binds
!>     the appropriate concrete implementation (Load or Excited).
!>   - `Initialize`: Block-assembly routine that iterates over body-type pairs
!>     and calls `TwoRdmInit_InitializeBlock` for each, then packs into state.
!>
!> @par Assembly Strategy
!> The full 2-RDM is built by iterating over body-type combinations:
!> @code
!>   for bt1 in {1,2}:
!>     for bt2 in {1,2}:
!>       rdm2Block = InitializeBlock(bt1, bt2)
!>       rdm2[start1:end1, start2:end2, start1:end1, start2:end2] = rdm2Block
!> @endcode
!> The assembled tensor is then packed via `TwoRdm_FillTwoRdmStateFromRdm2`.
submodule(M_TwoRdmInit) S_TwoRdmInit

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Factory dispatcher: reads JSON and binds 2-RDM initializer.
!>
!> @details
!> Checks for the existence of `"twoRdmInit.load"` or `"twoRdmInit.excited"`
!> in the active JSON configuration and calls the corresponding Fabricate
!> routine. The `TwoRdmInit_Initialize` pointer is bound to the local
!> `Initialize` routine which handles block assembly.
!>
!> @par Priority Order
!>   1. `"twoRdmInit.load"` - Load from files
!>   2. `"twoRdmInit.excited"` - Build from 1-RDM products
!>
!> @note Will `error stop` if neither option is specified.
  module subroutine TwoRdmInit_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_TwoRdmInit_Load
    use M_TwoRdmInit_Excited

    call Say_Fabricate("twoRdmInit")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    TwoRdmInit_Initialize => Initialize

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("twoRdmInit.load")) then
      call TwoRdmInit_Load_Fabricate

    else if (Json_GetExistence("twoRdmInit.excited")) then
      call TwoRdmInit_Excited_Fabricate

    else
      error stop "twoRdmInit is missing one of: load, excited"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Assemble full 2-RDM by iterating over body-type blocks.
!>
!> @details
!> This routine implements a block-by-block assembly strategy:
!>   1. Allocates a full (nO x nO x nO x nO) tensor where nO = total orbitals
!>   2. Loops over all (bt1, bt2) pairs in {1,2} x {1,2}
!>   3. For each pair, allocates appropriately-sized block and calls
!>      `TwoRdmInit_InitializeBlock` (bound by the concrete implementation)
!>   4. Copies block into the appropriate slice of the full tensor
!>   5. Packs the assembled tensor into the linearized state format
!>
!> @par Memory Considerations
!> The full tensor requires O(nOrbs^4) storage. For large orbital spaces,
!> consider using direct packing strategies in the TwoRdm module variants.
!>
!> @param[out] twoRdmState  Linearized 2-RDM state vector
  subroutine Initialize(twoRdmState)
    use M_Utils_Say
    use M_Utils_DataStorage
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_TwoRdm

    complex(R64), intent(out), contiguous  :: twoRdmState(:)

    complex(R64), allocatable :: rdm2(:, :, :, :)
    complex(R64), allocatable :: rdm2Block(:, :, :, :)
    integer(I32) :: nO
    integer(I32) :: ibt1, ibt2, startBt1, endBt1, startBt2, endBt2

    nO = Method_Mb_OrbBased_nOrbsSum

    allocate (rdm2(nO, nO, nO, nO))

    do ibt1 = 1, 2
      do ibt2 = 1, 2

        startBt1 = Method_Mb_OrbBased_nOrbsStart(ibt1)
        endBt1 = Method_Mb_OrbBased_nOrbsEnd(ibt1)
        startBt2 = Method_Mb_OrbBased_nOrbsStart(ibt2)
        endBt2 = Method_Mb_OrbBased_nOrbsEnd(ibt2)

        allocate (rdm2Block(endBt1 - startBt1 + 1, &
                            endBt2 - startBt2 + 1, &
                            endBt1 - startBt1 + 1, &
                            endBt2 - startBt2 + 1))

        call TwoRdmInit_InitializeBlock(rdm2Block, ibt1, ibt2)

        rdm2(startBt1:endBt1, startBt2:endBt2, &
             startBt1:endBt1, startBt2:endBt2) = rdm2Block

        deallocate (rdm2Block)

      end do
    end do

    call TwoRdm_FillTwoRdmStateFromRdm2(twoRdmState, rdm2)

  end subroutine

end submodule
