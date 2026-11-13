! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Singlet Valdemoro closure.
!>
!> @details
!> The closure Reconstruction_Singlet_FillRdm3aab: the free-cumulant part of
!> rdm3aab (3-cumulant zero). Driver: S_Reconstruction_Singlet.
submodule(M_Reconstruction_Singlet_Valdemoro) S_Reconstruction_Singlet_Valdemoro

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Singlet_Valdemoro_Fabricate
    use M_Utils_Say
    use M_Reconstruction_Singlet

    call Say_Fabricate("reconstruction.singlet.valdemoro")

    Reconstruction_Singlet_FillRdm3aab => FillRdm3aab

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Free-cumulant part (3-cumulant set to zero).
  subroutine FillRdm3aab(rdm3aab, rdm1a, rdm2ab)
    use M_Utils_RdmCumulants

    complex(R64), intent(out), allocatable :: rdm3aab(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm1a(:, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    call RdmCumulants_FillCumulFreePartM3aab(rdm3aab, rdm1a, rdm2ab)

  end subroutine

end submodule
