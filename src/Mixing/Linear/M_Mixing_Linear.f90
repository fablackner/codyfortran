! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Linear (damped) mixing backend.
!>
!> @details
!> Implements the classic damped fixed-point update
!>
!>   x ← (1−α)·x + α·xNew
!>
!> Stateless: no history is kept.
!>
!> **JSON Configuration:**
!> @code{.json}
!> "mixing": { "linear": { "alpha": 1.0 } }
!> @endcode
!>
!> | Parameter | Type | Default | Description                    |
!> |-----------|------|---------|--------------------------------|
!> | `alpha`   | real | 1.0     | Damping parameter (0 < α ≤ 1) |
!>
!> This backend is also selected when the JSON has no `mixing` block at all.
module M_Mixing_Linear
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface Mixing_Linear_Fabricate
    !> @brief Bind the linear mixing implementation to the Mixing pointer.
    module subroutine Mixing_Linear_Fabricate()
    end subroutine
  end interface

end module
