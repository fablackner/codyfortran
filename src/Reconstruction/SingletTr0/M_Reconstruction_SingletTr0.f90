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
!> suitable for constrained propagation and certain ML approaches.
!>
!> **Available schemes**:
!> - `valdemoro`: Free-cumulant approximation on trace-free representation
!> - `ny`: Nakatsuji-Yasuda correction on trace-free representation
!> - `mlpFullCop`: Neural network on packed trace-free features
!> - `cnnFullCop`: Convolutional neural network on 4D tensor input
!> - `cnnCorrelCop`: CNN correction to Valdemoro collision operator
!>
!> @see M_TwoRdm_Fermionic_SingletTr0 for state vector layout
module M_Reconstruction_SingletTr0
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Attach the selected SingletTr0 reconstruction implementation.
    !>
    !> Called during initialization to wire the procedure pointers exported by
    !> this module to a specific SingletTr0-aware scheme and to perform any
    !> lightweight common setup. Scheme-specific heavy setup should be done in
    !> the scheme’s own setup procedure if needed.
    module subroutine Reconstruction_SingletTr0_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
