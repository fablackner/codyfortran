! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of Singlet Valdemoro reconstruction.
!>
!> @details
!> **Algorithm**:
!> 1. Unpack rdm2ab from state vector
!> 2. Compute rdm1a via partial trace
!> 3. Reconstruct rdm3aab using cumulant-free expansion
!> 4. (Optional) Apply contraction-consistent trace constraints if ccQ=true
!> 5. Compute collision operator collOp2ab = Tr₃[h2, rdm3aab]
!> 6. Add to derivative and repack
submodule(M_Reconstruction_Singlet_Valdemoro) S_Reconstruction_Singlet_Valdemoro

  implicit none

  !> Flag: apply contraction-consistent (cc) trace constraints to rdm3aab
  logical :: ccQ

  !=============================================================================
  ! local data
  !=============================================================================

  !> Reconstructed 3-RDM α-α-β block; persists across calls
  complex(R64), allocatable :: rdm3aab(:, :, :, :, :, :)
  !> Collision operator contribution to 2-RDM α-β block derivative
  complex(R64), allocatable :: collOp2ab(:, :, :, :)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Singlet_Valdemoro_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Reconstruction

    call Say_Fabricate("reconstruction.singlet.singletValdemoro")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ccQ = Json_Get("reconstruction.singlet.valdemoro.ccQ", .false.)

    Reconstruction_AddCollisionOp => AddCollisionOp

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AddCollisionOp(dTwoRdmState, twoRdmState, h2)
    use M_Utils_RdmCumulants
    use M_Utils_RdmCollisionOp
    use M_Utils_RdmSingletRelations
    use M_Utils_RdmTrace
    use M_Utils_RdmEnforceTraceM3aab
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_TwoRdm_Fermionic_Singlet

    complex(R64), intent(inout), contiguous  :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)
    complex(R64), allocatable :: rdm1a(:, :)
    complex(R64), allocatable :: dRdm2ab(:, :, :, :)
    complex(R64), allocatable :: onefoldTraces(:, :, :, :, :)
    complex(R64), allocatable :: rdm2aa(:, :, :, :)

    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(dRdm2ab, dTwoRdmState)

    call RdmTrace_FillRdm1aFromRdm2ab(rdm1a, rdm2ab)

    call RdmCumulants_FillCumulFreePartM3aab(rdm3aab, rdm1a, rdm2ab)

    if (ccQ) then

      allocate (onefoldTraces(Method_Mb_OrbBased_nOrbs(1), &
                              Method_Mb_OrbBased_nOrbs(1), &
                              Method_Mb_OrbBased_nOrbs(1), &
                              Method_Mb_OrbBased_nOrbs(1), &
                              4))

      call RdmSingletRelations_FillRdm2aaFromRdm2ab(rdm2aa, rdm2ab)

      onefoldTraces(:, :, :, :, 1) = (Method_Mb_nBodies(1) - 1) * rdm2ab(:, :, :, :)
      onefoldTraces(:, :, :, :, 2) = rdm2aa(:, :, :, :)
      onefoldTraces(:, :, :, :, 3) = rdm2aa(:, :, :, :)
      onefoldTraces(:, :, :, :, 4) = (Method_Mb_nBodies(1)) * rdm2aa(:, :, :, :)

      call RdmEnforceTraceM3aab_Apply(rdm3aab, onefoldTraces)

    end if

    call RdmCollisionOp_FillCollisionOp2ab(collOp2ab, h2, rdm3aab)

    dRdm2ab(:, :, :, :) = dRdm2ab(:, :, :, :) + collOp2ab(:, :, :, :)

    call TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab(dTwoRdmState, dRdm2ab)

  end subroutine

end submodule
