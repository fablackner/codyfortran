! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Backend selector and wiring for fermionic 2‑RDM implementations.
!>
!> @details
!> This module serves as the dispatch layer for fermionic 2-RDM backends.
!> It does not implement algebra itself. Instead, it inspects the JSON
!> configuration and activates the appropriate symmetry-specialized backend.
!>
!> ## Available Fermionic Backends
!>
!> | Backend      | JSON Key                       | Storage Size              | Use Case                        |
!> |--------------|--------------------------------|---------------------------|---------------------------------|
!> | Generic      | `twoRdm.fermionic.generic`     | nOrbs⁴                    | Unrestricted, no spin symmetry  |
!> | Singlet      | `twoRdm.fermionic.singlet`     | nα² × nβ²                 | Spin-singlet, restricted        |
!> | SingletTr0   | `twoRdm.fermionic.singletTr0`  | nα² + nα² × nβ²           | Singlet with explicit 1-RDM     |
!>
!> The Normalize procedure is shared across all fermionic backends and
!> enforces the particle-number normalization: Tr(D²) = N(N-1).
!>
!> Use `TwoRdm_Fermionic_Fabricate` as part of the top‑level `TwoRdm_Fabricate`
!> to activate a fermionic backend.
module M_TwoRdm_Fermionic
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Select and activate a fermionic 2‑RDM backend.
    !>
    !> Responsibilities:
    !> - Inspect configuration and choose between available fermionic backends
    !>   (e.g. generic, singlet).
    !> - Assign all procedure pointers exported by `M_TwoRdm` to the chosen backend.
    !> - Optionally call the backend's setup to prepare storage/indexing.
    module subroutine TwoRdm_Fermionic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! local data
  !=============================================================================

end module

