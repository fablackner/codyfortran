! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Generic CI coefficient backend (representation-agnostic).
!>
!> This module provides the runtime wiring for a generic coefficient
!> representation that does not assume a particular model (e.g., Hubbard).
!> It connects the abstract interface in `M_Coeffs` to concrete, model-neutral
!> implementations based on the active configuration/basis description and
!> available body types (particle species, statistics).
module M_Coeffs_Generic
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Attach generic implementations to `M_Coeffs` and size the coefficient space.
    !>
    !> Responsibilities:
    !> - determine `Coeffs_nCoeffs` from the active configuration description,
    !> - associate procedure pointers in `M_Coeffs` with generic implementations,
    !> - perform any light-weight setup that is independent of a specific model.
    !>
    !> Heavy-weight, model-specific data (e.g., integrals) are not owned here and
    !> must be provided via arguments of the abstract interfaces when needed.
    module subroutine Coeffs_Generic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  ! (all exported via `M_Coeffs` after fabrication)
  !=============================================================================

end module

