! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of SingletTr0 Valdemoro reconstruction.
!>
!> @details
!> **Algorithm**:
!> 1. Unpack rdm2ab from trace-free state vector
!> 2. Compute rdm1a via partial trace
!> 3. Reconstruct rdm3aab using free-cumulant expansion
!> 4. Compute collOp2ab = Tr₃[h2, rdm3aab] and collOp1a = Tr₂[h2, rdm2ab]
!> 5. Project collOp2ab to trace-free subspace via RdmEnforceTraceM2ab_Apply
!> 6. Pack both contributions into state derivative vector
submodule(M_Reconstruction_SingletTr0_Valdemoro) S_Reconstruction_SingletTr0_Valdemoro

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Reconstructed 3-RDM α-α-β block; persists across calls
  complex(R64), allocatable :: rdm3aab(:, :, :, :, :, :)
  !> Collision operator contribution to 2-RDM α-β block derivative
  complex(R64), allocatable :: collOp2ab(:, :, :, :)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_SingletTr0_Valdemoro_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Reconstruction

    call Say_Fabricate("reconstruction.singletTr0.singletTr0Valdemoro")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Reconstruction_AddCollisionOp => AddCollisionOp

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AddCollisionOp(dTwoRdmState, twoRdmState, h2)
    use M_Utils_RdmCumulants
    use M_Utils_RdmEnforceTraceM2ab
    use M_Utils_RdmTrace
    use M_Utils_RdmCollisionOp
    use M_Method_Mb_OrbBased
    use M_TwoRdm_Fermionic_SingletTr0

    complex(R64), intent(inout), contiguous  :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)
    complex(R64), allocatable :: rdm1a(:, :)

    complex(R64), allocatable :: collOp1a(:, :)
    complex(R64), allocatable :: updateVec(:)

    integer(I32) :: nOa, nOb

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)

    call TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)

    call RdmTrace_FillRdm1aFromRdm2ab(rdm1a, rdm2ab)

    call RdmCumulants_FillCumulFreePartM3aab(rdm3aab, rdm1a, rdm2ab)

    call RdmCollisionOp_FillCollisionOp2ab(collOp2ab, h2, rdm3aab)
    call RdmCollisionOp_FillCollisionOp1a(collOp1a, h2, rdm2ab)

    call RdmEnforceTraceM2ab_Apply(collOp2ab)

    allocate (updateVec(size(dTwoRdmState)))
    call TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdms(updateVec, rdm2abTr0_=collOp2ab, rdm1a_=collOp1a)

    dTwoRdmState = dTwoRdmState + updateVec

  end subroutine

end submodule
