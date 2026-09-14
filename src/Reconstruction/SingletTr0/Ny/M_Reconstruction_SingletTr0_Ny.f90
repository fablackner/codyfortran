! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief SingletTr0 Nakatsuji-Yasuda (NY) reconstruction scheme.
!>
!> @details
!> Extends the Valdemoro approximation by adding the NY correction to the
!> 3-body cumulant, operating on the trace-free singlet state representation.
!>
!> **JSON path**: `reconstruction.singletTr0.ny`
!>
!> @see H. Nakatsuji & K. Yasuda, Phys. Rev. Lett. 76, 1039 (1996)
module M_Reconstruction_SingletTr0_Ny
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the singletTr0 NY implementation at runtime.
    module subroutine Reconstruction_SingletTr0_Ny_Fabricate
    end subroutine
  end interface

end module
