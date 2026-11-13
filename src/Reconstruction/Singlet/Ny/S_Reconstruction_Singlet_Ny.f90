! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Singlet Nakatsuji-Yasuda closure.
!>
!> @details
!> The closure Reconstruction_Singlet_FillRdm3aab: the free-cumulant part plus
!> the NY 3-cumulant, a sum of products of two 2-cumulants contracted through
!> R = 2 P_occ - 1 (P_occ projects on the nBodies(1) most strongly occupied
!> natural orbitals of rdm1a; R depends on the natural orbitals only). Driver:
!> S_Reconstruction_Singlet.
submodule(M_Reconstruction_Singlet_Ny) S_Reconstruction_Singlet_Ny

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Singlet_Ny_Fabricate
    use M_Utils_Say
    use M_Reconstruction_Singlet

    call Say_Fabricate("reconstruction.singlet.singletNy")

    Reconstruction_Singlet_FillRdm3aab => FillRdm3aab

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Free-cumulant part plus the NY 3-cumulant.
  subroutine FillRdm3aab(rdm3aab, rdm1a, rdm2ab)
    use M_Utils_RdmCumulants
    use M_Method_Mb

    complex(R64), intent(out), allocatable :: rdm3aab(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm1a(:, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    call RdmCumulants_FillCumulFreePartM3aab(rdm3aab, rdm1a, rdm2ab)
    call RdmCumulants_AddNyCumulantM3aab(rdm3aab, Method_Mb_nBodies(1), rdm1a, rdm2ab)

  end subroutine

end submodule
