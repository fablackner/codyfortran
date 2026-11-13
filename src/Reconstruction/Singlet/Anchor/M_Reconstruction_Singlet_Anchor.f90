! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Anchored singlet closure (`reconstruction.singlet.anchor`).
!>
!> @details
!> An extension of the singlet family: it wraps the closure of the scheme by
!> the closure error at a reference state, and with `nTransitions` > 0 the
!> collision operator by the spectral part. The method, the family driver and
!> the schemes do not know about it. Details:
!> src/Reconstruction/Singlet/Anchor/A_Reconstruction_Singlet_Anchor.md.
module M_Reconstruction_Singlet_Anchor
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Reads reconstruction.singlet.anchor and wraps the closure (and the
    !> collision operator) bound by the singlet family; binds
    !> Reconstruction_Setup and Reconstruction_SetOrbitals.
    module subroutine Reconstruction_Singlet_Anchor_Fabricate
    end subroutine

    !> Writes the reference files of the anchor from CI vectors, ket-first
    !> with the transition index last: rdm3Bt01_01_02.in, the up-up-down
    !> 3-RDM of kets(:, 1), and with more than one ket trdm3Bt01_01_02.in and
    !> trdm2Bt01_02.in, the transition RDMs <kets(:,1)|..|kets(:,k)>, k >= 2
    !> (the reference 2-RDM rdm2Bt01_02.in comes from Coeffs_SaveTwoRdm).
    module subroutine Reconstruction_Singlet_Anchor_SaveReference(kets)
      complex(R64), intent(in), contiguous :: kets(:, :)
    end subroutine
  end interface

end module
