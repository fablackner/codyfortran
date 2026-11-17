! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Singlet/Reduce purification backend.
!>
!> @details
!> Implements iterative positivity reduction for systems with singlet spin symmetry.
!> The algorithm exploits the spin-adapted block structure of the 2-RDM to separately
!> purify singlet (symmetric) and triplet (antisymmetric) components.
!>
!> @par Mathematical Foundation
!> For a singlet state, the 2-RDM in the α-β sector decomposes as:
!> @f[
!>   D^2_{\alpha\beta} = D^2_{\text{singlet}} + D^2_{\text{triplet}}
!> @f]
!> where singlet (S=0) couples symmetrically and triplet (S=1) antisymmetrically
!> under particle exchange. The same decomposition applies to the 2-hole-RDM Q².
!>
!> @par Iterative Reduction Algorithm
!> @verbatim
!>   1. Unpack rdm2ab from twoRdmState
!>   2. LOOP until min(eigenvalues) > threshold OR iteration > nSteps:
!>        a. Split rdm2ab → singlet + triplet blocks
!>        b. Eigendecompose each block; extract defective (λ < threshold) part
!>        c. Make defects trace-free (contraction-free correction)
!>        d. Subtract defects from rdm2ab
!>        e. Compute qrdm2ab = Q-RDM from updated rdm2ab
!>        f. Repeat (a-d) for Q-RDM blocks
!>        g. Optionally project onto energy-conserving subspace (⟂ h²)
!>   3. Repack rdm2ab → twoRdmState
!> @endverbatim
!>
!> @par Configuration Parameters
!> | JSON Path                                    | Default   | Description                              |
!> |----------------------------------------------|-----------|------------------------------------------|
!> | `purification.singlet.reduce.nSteps`         | 40        | Maximum purification iterations          |
!> | `purification.singlet.reduce.threshold`      | -1e-5     | Eigenvalue threshold (negative allowed)  |
!> | `purification.singlet.reduce.conserveEnergyQ`| .FALSE.   | Project corrections ⊥ to h² (energy cons)|
!>
!> @par Effects
!> Binds `M_Purification::Purification_PurifyTwoRdmState` to the local implementation.
!>
!> @see M_TwoRdm_Fermionic_Singlet for the packed state ↔ rdm2ab conversion
!> @see M_Utils_RdmDefective for the defective eigenspace extraction
!> @see M_Utils_RdmHole for Q-RDM construction
module M_Purification_Singlet_Reduce
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Read parameters for the reduce backend and assign the global
    !> `M_Purification` procedure pointer `Purification_PurifyTwoRdmState`.
    module subroutine Purification_Singlet_Reduce_Fabricate
    end subroutine

    !> Trace-free direction of the two-body energy of a singlet as a
    !> functional of its (up,down) 2-RDM block D_ab: with D_aa = D_ab - P12 D_ab
    !> the energy is <2 h2_singlet + 6 h2_triplet | D_ab>, so corrections
    !> orthogonal to h2Tr0 = Tr0(h2_singlet) + 3 Tr0(h2_triplet) of the
    !> (up,down) block of h2 conserve it. Shared with SingletTr0/Reduce.
    module subroutine Purification_Singlet_Reduce_FillEnergyDirection(h2Tr0, h2, nOa)
      !> Output: the direction, (nOa, nOa, nOa, nOa).
      complex(R64), intent(out), allocatable :: h2Tr0(:, :, :, :)
      !> Two-body Hamiltonian, the full spin-orbital tensor (2 nOa)^4.
      complex(R64), intent(in), contiguous :: h2(:, :, :, :)
      !> Number of spin-up orbitals.
      integer(I32), intent(in) :: nOa
    end subroutine

    !> The iterative positivity reduction of a singlet (up,down) 2-RDM block,
    !> shared with SingletTr0/Reduce. Every step removes the defective parts
    !> (eigenvalues below threshold) of the singlet and triplet blocks of D,
    !> then those of Q of the updated D; the corrections are trace-free and,
    !> with h2Tr0_, orthogonal to the energy direction
    !> (Purification_Singlet_Reduce_FillEnergyDirection). Stops when the
    !> smallest eigenvalue reaches threshold or after nSteps steps.
    module subroutine Purification_Singlet_Reduce_PurifyRdm2ab(rdm2ab, nSteps, threshold, h2Tr0_)
      !> The (up,down) 2-RDM block, purified in place.
      complex(R64), intent(inout), contiguous :: rdm2ab(:, :, :, :)
      !> Largest number of steps.
      integer(I32), intent(in) :: nSteps
      !> Eigenvalue threshold of the defect detection.
      real(R64), intent(in) :: threshold
      !> Energy direction of energy-conserving corrections.
      complex(R64), intent(in), contiguous, optional :: h2Tr0_(:, :, :, :)
    end subroutine
  end interface

  !=============================================================================
  ! local data
  !=============================================================================

end module
