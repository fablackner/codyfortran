! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief SingletTr0/Reduce purification backend.
!>
!> @details
!> Implements iterative positivity reduction for systems using the SingletTr0
!> 2-RDM representation, which stores only the trace-free (cumulant-like)
!> component of the α-β 2-RDM block alongside an explicit 1-RDM.
!>
!> @par Key Difference from Singlet Backend
!> In the SingletTr0 representation, the full 2-RDM is reconstructed as:
!> @f[
!>   D^2_{\alpha\beta}(pq|rs) = D^{2,\text{Tr0}}_{\alpha\beta} + D^1_\alpha(pr) D^1_\beta(qs)
!> @f]
!> This reduces the DOF count but requires the trace constraint to be enforced
!> via `RdmEnforceTraceM2ab` rather than the full symmetric/antisymmetric projectors.
!>
!> @par Configuration Parameters
!> | JSON Path                                       | Default   | Description                              |
!> |-------------------------------------------------|-----------|------------------------------------------|
!> | `purification.singletTr0.reduce.nSteps`         | 40        | Maximum purification iterations          |
!> | `purification.singletTr0.reduce.threshold`      | -1e-5     | Eigenvalue threshold (negative allowed)  |
!> | `purification.singletTr0.reduce.conserveEnergyQ`| .FALSE.   | Project corrections ⊥ to h² (energy cons)|
!>
!> @par Effects
!> Binds `M_Purification::Purification_PurifyTwoRdmState` to the local implementation.
!>
!> @see M_TwoRdm_Fermionic_SingletTr0 for the packed state ↔ rdm2ab conversion
!> @see M_Purification_Singlet_Reduce for the analogous full-trace algorithm
module M_Purification_SingletTr0_Reduce
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Read parameters for the reduce backend and assign the global
    !> `M_Purification` procedure pointer `Purification_PurifyTwoRdmState`.
    module subroutine Purification_SingletTr0_Reduce_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! local data
  !=============================================================================

end module
