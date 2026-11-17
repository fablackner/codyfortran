! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of Generic Valdemoro reconstruction.
!>
!> @details
!> **Algorithm**:
!> 1. Unpack 2-RDM and 1-RDM from state vector
!> 2. Reconstruct 3-RDM using `RdmCumulants_FillCumulFreePartM3` (C³=0 approximation)
!> 3. Compute collision operator collOp2 = Tr₃[h2, rdm3]
!> 4. Add collision term to derivative and repack
submodule(M_Reconstruction_Generic_Valdemoro) S_Reconstruction_Generic_Valdemoro

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Generic_Valdemoro_Fabricate
    use M_Utils_Say
    use M_Reconstruction

    call Say_Fabricate("reconstruction.generic.valdemoro")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Reconstruction_AddCollisionOp => AddCollisionOp

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AddCollisionOp(dTwoRdmState, twoRdmState, h2)
    use M_Utils_RdmCumulants
    use M_TwoRdm
    use M_Utils_RdmCollisionOp

    complex(R64), intent(inout), contiguous  :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)

    complex(R64), allocatable :: rdm2(:, :, :, :)
    complex(R64), allocatable :: rdm1(:, :)
    complex(R64), allocatable :: dRdm2(:, :, :, :)
    complex(R64), allocatable :: rdm3(:, :, :, :, :, :)
    complex(R64), allocatable :: collOp2(:, :, :, :)

    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, twoRdmState)
    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)
    call TwoRdm_FillRdm2FromTwoRdmState(dRdm2, dTwoRdmState)

    call RdmCumulants_FillCumulFreePartM3(rdm3, rdm1, rdm2)

    call RdmCollisionOp_FillCollisionOp2(collOp2, h2, rdm3)

    dRdm2(:, :, :, :) = dRdm2(:, :, :, :) + collOp2(:, :, :, :)

    call TwoRdm_FillTwoRdmStateFromRdm2(dTwoRdmState, dRdm2)

  end subroutine

end submodule
