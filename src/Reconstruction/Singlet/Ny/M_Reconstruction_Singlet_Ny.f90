! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Singlet Nakatsuji-Yasuda (NY) reconstruction scheme.
!>
!> @details
!> Extends the Valdemoro approximation by the NY 3-cumulant: products of two
!> 2-cumulants contracted through R = 2 P_occ - 1.
!>
!> **JSON path**: `reconstruction.singlet.singletNy`
!>
!> **Options** (read by the family, S_Reconstruction_Singlet):
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
    !> Binds the closure Reconstruction_Singlet_FillRdm3aab.
    module subroutine Reconstruction_Singlet_Ny_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! local data
  !=============================================================================

end module
