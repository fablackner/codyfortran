! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief SingletTr0 Valdemoro closure: the cumulant-free
!> part of rdm3aab, 3-cumulant zero (the driver is S_Reconstruction_SingletTr0).
submodule(M_Reconstruction_SingletTr0_Valdemoro) S_Reconstruction_SingletTr0_Valdemoro

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_SingletTr0_Valdemoro_Fabricate
    use M_Utils_Say
    use M_Reconstruction_SingletTr0

    call Say_Fabricate("reconstruction.singletTr0.valdemoro")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Reconstruction_SingletTr0_FillRdm3aab => FillRdm3aab

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm3aab(rdm3aab, rdm1a, rdm2ab)
    use M_Utils_RdmCumulants

    complex(R64), intent(out), allocatable :: rdm3aab(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm1a(:, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    call RdmCumulants_FillCumulFreePartM3aab(rdm3aab, rdm1a, rdm2ab)

  end subroutine

end submodule
