! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Singlet family: dispatch to the scheme and the common driver; see
!> M_Reconstruction_Singlet.
submodule(M_Reconstruction_Singlet) S_Reconstruction_Singlet

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Singlet_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Reconstruction
    use M_Reconstruction_Singlet_Ny
    use M_Reconstruction_Singlet_Valdemoro
    use M_Reconstruction_Singlet_Triples
    use M_Reconstruction_Singlet_Anchor

    character(len=:), allocatable :: scheme

    call Say_Fabricate("reconstruction.singlet")

    Reconstruction_Singlet_FillRdm3aab => null()
    if (Json_GetExistence("reconstruction.singlet.singletNy")) then
      scheme = "reconstruction.singlet.singletNy"
      call Reconstruction_Singlet_Ny_Fabricate
    else if (Json_GetExistence("reconstruction.singlet.valdemoro")) then
      scheme = "reconstruction.singlet.valdemoro"
      call Reconstruction_Singlet_Valdemoro_Fabricate
    else
      error stop "reconstruction.singlet is missing one of: singletNy, valdemoro"
    end if
    if (.not. associated(Reconstruction_Singlet_FillRdm3aab)) error stop "reconstruction.singlet: the scheme binds no closure"
    Reconstruction_Singlet_ccQ = Json_Get(scheme//".ccQ", .false.)

    Reconstruction_AddCollisionOp => AddCollisionOp

    ! optional extensions, each wrapping the closure bound before it
    ! (Triples/, Anchor/)
    if (Json_GetExistence("reconstruction.singlet.triples")) call Reconstruction_Singlet_Triples_Fabricate
    if (Json_GetExistence("reconstruction.singlet.anchor")) call Reconstruction_Singlet_Anchor_Fabricate

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AddCollisionOp(dTwoRdmState, twoRdmState, h2)
    use M_Utils_RdmCollisionOp
    use M_TwoRdm_Fermionic_Singlet

    complex(R64), intent(inout), contiguous :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous    :: twoRdmState(:)
    complex(R64), intent(in), contiguous    :: h2(:, :, :, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :), dRdm2ab(:, :, :, :), rdm3aab(:, :, :, :, :, :)
    complex(R64), allocatable :: collOp2ab(:, :, :, :)

    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(dRdm2ab, dTwoRdmState)

    call Reconstruction_Singlet_ReconstructRdm3aab(rdm3aab, rdm2ab)
    call RdmCollisionOp_FillCollisionOp2ab(collOp2ab, h2, rdm3aab)
    dRdm2ab(:, :, :, :) = dRdm2ab(:, :, :, :) + collOp2ab(:, :, :, :)

    call TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab(dTwoRdmState, dRdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Singlet_ReconstructRdm3aab(rdm3aab, rdm2ab)
    use M_Utils_RdmTrace

    complex(R64), intent(out), allocatable :: rdm3aab(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    complex(R64), allocatable :: rdm1a(:, :)

    call RdmTrace_FillRdm1aFromRdm2ab(rdm1a, rdm2ab)
    call Reconstruction_Singlet_FillRdm3aab(rdm3aab, rdm1a, rdm2ab)
    if (Reconstruction_Singlet_ccQ) call EnforceContractions(rdm3aab, rdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Contraction consistency of rdm3aab with rdm2ab (ccQ).
  subroutine EnforceContractions(rdm3aab, rdm2ab)
    use M_Utils_RdmSingletRelations
    use M_Utils_RdmEnforceTraceM3aab
    use M_Method_Mb

    complex(R64), intent(inout), contiguous :: rdm3aab(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous    :: rdm2ab(:, :, :, :)

    complex(R64), allocatable :: onefoldTraces(:, :, :, :, :), rdm2aa(:, :, :, :)
    integer(I32) :: n, nUp

    n = size(rdm2ab, 1)
    nUp = Method_Mb_nBodies(1)

    allocate (onefoldTraces(n, n, n, n, 4))
    call RdmSingletRelations_FillRdm2aaFromRdm2ab(rdm2aa, rdm2ab)

    onefoldTraces(:, :, :, :, 1) = (nUp - 1) * rdm2ab(:, :, :, :)
    onefoldTraces(:, :, :, :, 2) = rdm2aa(:, :, :, :)
    onefoldTraces(:, :, :, :, 3) = rdm2aa(:, :, :, :)
    onefoldTraces(:, :, :, :, 4) = nUp * rdm2aa(:, :, :, :)

    call RdmEnforceTraceM3aab_Apply(rdm3aab, onefoldTraces_=onefoldTraces)

  end subroutine

end submodule
