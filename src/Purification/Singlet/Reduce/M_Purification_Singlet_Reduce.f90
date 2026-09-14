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
  end interface

  !=============================================================================
  ! local data
  !=============================================================================

end module
