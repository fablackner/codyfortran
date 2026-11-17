! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Trace-free spin-singlet (SingletTr0) reconstruction contracts.
!>
!> @details
!> This module declares the interfaces for reconstruction schemes that operate
!> on **trace-free singlet state vectors**. In the SingletTr0 representation,
!> the 2-RDM is split into:
!>
!> - `rdm2abTr0`: Trace-free α-β block (reduced DOF)
!> - `rdm1a`: 1-RDM stored explicitly (provides the trace)
!>
!> This representation has fewer degrees of freedom than full singlet, making it
!> suitable for constrained propagation.
!>
!> **Available schemes**:
!> - `valdemoro`: Free-cumulant approximation on trace-free representation
!> - `ny`: Nakatsuji-Yasuda correction on trace-free representation
!>
!> A scheme fabricate binds only the closure
!> Reconstruction_SingletTr0_FillRdm3aab (the closures of the Singlet family,
!> I_Reconstruction_Singlet_FillRdm3aab); the family driver
!> (S_Reconstruction_SingletTr0) unpacks the state, adds the collision
!> operator of the 2-RDM and of the 1-RDM and packs the trace-free result.
!>
!> @see M_TwoRdm_Fermionic_SingletTr0 for state vector layout
module M_Reconstruction_SingletTr0
  use M_Utils_Types
  use M_Reconstruction_Singlet, only: I_Reconstruction_Singlet_FillRdm3aab

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds the closure of the scheme and the family driver
    !> (Reconstruction_AddCollisionOp).
    module subroutine Reconstruction_SingletTr0_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> The closure of the selected scheme: rdm3aab from rdm1a and a singlet
  !> rdm2ab.
  procedure(I_Reconstruction_Singlet_FillRdm3aab), pointer :: Reconstruction_SingletTr0_FillRdm3aab => null()

end module
