! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for M_TwoRdm_Fermionic.
!>
!> @details
!> Implements the fermionic backend fabrication and shared normalization logic.
!> Dispatches to one of: Generic, Singlet, or SingletTr0 based on JSON config.
submodule(M_TwoRdm_Fermionic) S_TwoRdm_Fermionic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Fabricate the fermionic 2-RDM backend from JSON configuration.
!>
!> @details
!> Reads the `twoRdm.fermionic` section and activates the appropriate
!> symmetry-specialized implementation. Sets the shared Normalize procedure.
  module subroutine TwoRdm_Fermionic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_TwoRdm
    use M_TwoRdm_Fermionic_Generic
    use M_TwoRdm_Fermionic_Singlet
    use M_TwoRdm_Fermionic_SingletTr0

    call Say_Fabricate("twoRdm.fermionic")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    TwoRdm_Normalize => Normalize

    !------------------------------------
    ! branch on spin symmetry
    !------------------------------------

    if (Json_GetExistence("twoRdm.fermionic.generic")) then
      call TwoRdm_Fermionic_Generic_Fabricate

    else if (Json_GetExistence("twoRdm.fermionic.singlet")) then
      call TwoRdm_Fermionic_Singlet_Fabricate

    else if (Json_GetExistence("twoRdm.fermionic.singletTr0")) then
      call TwoRdm_Fermionic_SingletTr0_Fabricate

    else
      error stop "twoRdm.fermionic is missing one of: generic, singlet, singletTr0"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Enforce particle-number normalization on the 2-RDM.
!>
!> @details
!> Scales the 2-RDM such that its trace equals N(N-1), where N is the total
!> particle number. This is a fundamental N-representability condition:
!>
!>   Tr(D²) = Σ_{i,j} D²(i,j;i,j) = N(N-1)
!>
!> The procedure unpacks the linearized state, computes the current trace,
!> rescales, and repacks. This is O(nOrbs⁴) due to the unpacking overhead
!> but ensures numerical stability during propagation.
!>
!> @param[inout] twoRdmState  Linearized 2-RDM to normalize in place
  subroutine Normalize(twoRdmState)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_TwoRdm

    complex(R64), intent(inout), contiguous  :: twoRdmState(:)

    complex(R64), allocatable :: rdm2(:, :, :, :)
    complex(R64) :: norm
    integer(I32) :: i1, i2, nO

    nO = Method_Mb_OrbBased_nOrbsSum

    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)

    ! Compute trace: Tr(D²) = Σ_{i,j} D²(i,j;i,j)
    norm = 0.0_R64
    do i1 = 1, nO
      do i2 = 1, nO

        norm = norm + rdm2(i1, i2, i1, i2)

      end do
    end do

    ! Scale to target normalization N(N-1)
    rdm2(:, :, :, :) = rdm2(:, :, :, :) / norm * (Method_Mb_nBodiesSum - 1) * Method_Mb_nBodiesSum

    call TwoRdm_FillTwoRdmStateFromRdm2(twoRdmState, rdm2)

  end subroutine

end submodule
