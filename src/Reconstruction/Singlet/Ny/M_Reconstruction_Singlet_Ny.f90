! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Singlet Nakatsuji-Yasuda (NY) reconstruction scheme.
!>
!> @details
!> Extends the Valdemoro approximation by adding particle-number-dependent
!> corrections to the 3-body cumulant. The NY correction improves accuracy
!> for systems where correlation is significant.
!>
!> **JSON path**: `reconstruction.singlet.singletNy`
!>
!> **Options**:
!> - `ccQ` (bool, default=false): Apply contraction-consistent trace constraints
!>
!> @see H. Nakatsuji & K. Yasuda, Phys. Rev. Lett. 76, 1039 (1996)
module M_Reconstruction_Singlet_Ny
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the singlet NY implementation at runtime.
    !>
    !> Wires the singlet-specific procedure pointers (e.g., alpha–beta block
    !> builders) to NY implementations and may read/apply variant options as
    !> required by the codebase’s configuration system.
    module subroutine Reconstruction_Singlet_Ny_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! local data
  !=============================================================================

end module
