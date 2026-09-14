! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Spin-unresolved (generic) reconstruction contracts.
!>
!> @details
!> This module declares the interfaces for generic reconstruction schemes that
!> operate on the **full 2-RDM tensor** without exploiting spin symmetry. The
!> RDM tensors have dimensions nOrbs⁴ (2-RDM) and nOrbs⁶ (3-RDM).
!>
!> **Available schemes**:
!> - `genericValdemoro`: Free-cumulant (Valdemoro) approximation with C³=0
!>
!> Concrete implementations live in submodules and are attached at runtime by
!> `Reconstruction_Generic_Fabricate`.
!>
!> @note Generic schemes are computationally more expensive than singlet variants
!>       due to the larger tensor dimensions. Use singlet schemes when applicable.
module M_Reconstruction_Generic
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Attach the selected generic reconstruction implementation.
    !>
    !> Called during initialization to wire the procedure pointers exported by
    !> this module to a specific scheme (e.g., Valdemoro). May also perform
    !> light-weight setup shared by all generic schemes.
    module subroutine Reconstruction_Generic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module

