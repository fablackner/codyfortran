! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Spin-singlet (S=0) reconstruction contracts.
!>
!> @details
!> For S=0 states only the α-β block rdm2ab(nα,nα,nβ,nβ) of the 2-RDM is
!> independent; the α-α and β-β blocks follow from it.
!>
!> **Schemes**: `valdemoro` (free-cumulant closure, 3-cumulant zero) and
!> `ny` (Nakatsuji-Yasuda 3-cumulant). A scheme fabricate binds only the
!> closure Reconstruction_Singlet_FillRdm3aab, never AddCollisionOp or Setup:
!> the family driver (S_Reconstruction_Singlet) unpacks the state, enforces
!> the contraction consistency (`<scheme>.ccQ`) and adds the collision
!> operator. Extensions wrap the closure and live in subfolders (Triples/,
!> Anchor/).
!>
!> All RDMs are ket-first (the layout of TwoRdm):
!>   rdm2ab(i1,i2,j1,j2)        = <a+_j1a a+_j2b a_i2b a_i1a>,
!>   rdm3aab(i1,i2,i3,j1,j2,j3) = <a+_j1a a+_j2a a+_j3b a_i3b a_i2a a_i1a>.
!>
!> @see M_TwoRdm_Fermionic_Singlet for state vector layout
module M_Reconstruction_Singlet
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds the closure of the scheme and the family driver.
    module subroutine Reconstruction_Singlet_Fabricate
    end subroutine

    !> The reconstructed 3-RDM of a singlet rdm2ab: the closure, then the
    !> contraction consistency if Reconstruction_Singlet_ccQ.
    module subroutine Reconstruction_Singlet_ReconstructRdm3aab(rdm3aab, rdm2ab)
      complex(R64), intent(out), allocatable :: rdm3aab(:, :, :, :, :, :)
      complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Contraction consistency of the reconstructed rdm3aab (`<scheme>.ccQ`,
  !> set by Reconstruction_Singlet_Fabricate).
  logical :: Reconstruction_Singlet_ccQ = .false.

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> The closure R of the selected scheme: every element of rdm3aab from rdm1a
  !> and a singlet rdm2ab (rdm2ab(i1,i2,j1,j2) = rdm2ab(i2,i1,j2,j1));
  !> antisymmetric in (i1,i2) and (j1,j2), not contraction-consistent.
  procedure(I_Reconstruction_Singlet_FillRdm3aab), pointer :: Reconstruction_Singlet_FillRdm3aab => null()
  abstract interface
    subroutine I_Reconstruction_Singlet_FillRdm3aab(rdm3aab, rdm1a, rdm2ab)
      import :: R64
      complex(R64), intent(out), allocatable :: rdm3aab(:, :, :, :, :, :)
      complex(R64), intent(in), contiguous   :: rdm1a(:, :)
      complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)
    end subroutine
  end interface

end module
