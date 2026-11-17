! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the SingletTr0/Reduce purification backend.
submodule(M_Purification_SingletTr0_Reduce) S_Purification_SingletTr0_Reduce

  implicit none

  !---------------------------------------------------------------------------
  ! Backend configuration (module-scope, set during fabrication)
  !---------------------------------------------------------------------------
  real(R64)    :: threshold       !< Eigenvalue cutoff for defect detection
  integer(I32) :: nSteps          !< Maximum iterations before forced exit
  logical      :: conserveEnergyQ !< If .TRUE., project corrections ⊥ to h²

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Configure the SingletTr0/Reduce backend and bind procedure pointers.
  module subroutine Purification_SingletTr0_Reduce_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Purification

    call Say_Fabricate("purification.singletTr0.reduce")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    nSteps = Json_Get("purification.singletTr0.reduce.nSteps", 40)
    threshold = Json_Get("purification.singletTr0.reduce.threshold", -1e-5_R64)
    conserveEnergyQ = Json_Get("purification.singletTr0.reduce.conserveEnergyQ", .FALSE.)
    Purification_usesH2Q = conserveEnergyQ

    Purification_PurifyTwoRdmState => PurifyTwoRdmState

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Purify a SingletTr0-symmetry 2-RDM state by iterative defect reduction.
  !>
  !> @details
  !> Unpacks the (up,down) block from the trace-free representation, reduces
  !> it with the kernel of Singlet/Reduce (Purification_Singlet_Reduce_PurifyRdm2ab:
  !> D corrections, then the Q corrections of the updated D) and packs it.
  !> The corrections are trace-free, so the 1-RDM segment is unchanged.
  !>
  !> @param[inout] twoRdmState  Packed 2-RDM (SingletTr0 layout).
  !> @param[in]    h2           Two-body Hamiltonian for optional energy
  !>                            conservation (full spin-orbital tensor).
  subroutine PurifyTwoRdmState(twoRdmState, h2)
    use M_TwoRdm_Fermionic_SingletTr0
    use M_Purification_Singlet_Reduce, only: Purification_Singlet_Reduce_FillEnergyDirection, &
                                             Purification_Singlet_Reduce_PurifyRdm2ab

    complex(R64), intent(inout), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous     :: h2(:, :, :, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)        !< Full α-β 2-RDM
    complex(R64), allocatable :: h2Tr0(:, :, :, :)         !< Trace-free energy direction

    call TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)

    if (conserveEnergyQ) then
      call Purification_Singlet_Reduce_FillEnergyDirection(h2Tr0, h2, size(rdm2ab, 1))
      call Purification_Singlet_Reduce_PurifyRdm2ab(rdm2ab, nSteps, threshold, h2Tr0_=h2Tr0)
    else
      call Purification_Singlet_Reduce_PurifyRdm2ab(rdm2ab, nSteps, threshold)
    end if

    call TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdm2ab(twoRdmState, rdm2ab)

  end subroutine

end submodule
