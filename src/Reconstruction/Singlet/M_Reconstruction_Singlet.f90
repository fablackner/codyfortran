! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Spin-singlet (S=0) reconstruction contracts.
!>
!> @details
!> This module declares the interfaces for reconstruction schemes that exploit
!> **spin-singlet symmetry**. For S=0 states, the 2-RDM has a block structure
!> where only the α-β block rdm2ab(nα,nα,nβ,nβ) is independent; α-α and β-β
!> blocks can be derived from it.
!>
!> **Available schemes**:
!> - `valdemoro`   : Free-cumulant (C³=0) approximation
!> - `singletNy`   : Nakatsuji-Yasuda with particle-number corrections
!> - `nnCholesky`  : Neural network using Cholesky-decomposed features
!>
!> The singlet structure reduces memory and compute by ~4× compared to generic schemes.
!>
!> @see M_TwoRdm_Fermionic_Singlet for state vector layout
module M_Reconstruction_Singlet
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Attach the selected singlet reconstruction implementation.
    !>
    !> Called during initialization to wire the procedure pointers exported by
    !> this module to a specific singlet-aware scheme and to perform any
    !> lightweight common setup. Scheme-specific heavy setup should be done in
    !> the scheme’s own setup procedure if needed.
    module subroutine Reconstruction_Singlet_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module

