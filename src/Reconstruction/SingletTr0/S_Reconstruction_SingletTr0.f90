! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> SingletTr0 family: dispatch to the scheme and the common driver; see
!> M_Reconstruction_SingletTr0.
submodule(M_Reconstruction_SingletTr0) S_Reconstruction_SingletTr0

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_SingletTr0_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Reconstruction
    use M_Reconstruction_SingletTr0_Valdemoro
    use M_Reconstruction_SingletTr0_Ny

    call Say_Fabricate("reconstruction.singletTr0")

    !------------------------------------
    ! branch
    !------------------------------------

    Reconstruction_SingletTr0_FillRdm3aab => null()

    if (Json_GetExistence("reconstruction.singletTr0.valdemoro")) then
      call Reconstruction_SingletTr0_Valdemoro_Fabricate

    else if (Json_GetExistence("reconstruction.singletTr0.ny")) then
      call Reconstruction_SingletTr0_Ny_Fabricate

    else
      error stop "reconstruction.singletTr0 is missing one of: valdemoro, ny"
    end if

    if (.not. associated(Reconstruction_SingletTr0_FillRdm3aab)) &
      error stop "reconstruction.singletTr0: the scheme binds no closure"
    Reconstruction_AddCollisionOp => AddCollisionOp

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> dTwoRdmState += the collision operators of the 2-RDM, Tr3[h2, rdm3aab]
  !> projected on the trace-free subspace, and of the 1-RDM, Tr2[h2, rdm2ab],
  !> with rdm3aab of the closure.
  subroutine AddCollisionOp(dTwoRdmState, twoRdmState, h2)
    use M_Utils_RdmEnforceTraceM2ab
    use M_Utils_RdmTrace
    use M_Utils_RdmCollisionOp
    use M_Method_Mb
    use M_TwoRdm_Fermionic_SingletTr0

    complex(R64), intent(inout), contiguous :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous    :: twoRdmState(:)
    complex(R64), intent(in), contiguous    :: h2(:, :, :, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :), rdm1a(:, :), rdm3aab(:, :, :, :, :, :)
    complex(R64), allocatable :: collOp2ab(:, :, :, :), collOp1a(:, :), updateVec(:)

    call TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    call RdmTrace_FillRdm1aFromRdm2ab(rdm1a, rdm2ab, Method_Mb_nBodies(2))

    call Reconstruction_SingletTr0_FillRdm3aab(rdm3aab, rdm1a, rdm2ab)

    call RdmCollisionOp_FillCollisionOp2ab(collOp2ab, h2, rdm3aab)
    call RdmCollisionOp_FillCollisionOp1a(collOp1a, h2, rdm2ab)

    call RdmEnforceTraceM2ab_Apply(collOp2ab)

    allocate (updateVec(size(dTwoRdmState)))
    call TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdms(updateVec, rdm2abTr0_=collOp2ab, rdm1a_=collOp1a)

    dTwoRdmState = dTwoRdmState + updateVec

  end subroutine

end submodule
