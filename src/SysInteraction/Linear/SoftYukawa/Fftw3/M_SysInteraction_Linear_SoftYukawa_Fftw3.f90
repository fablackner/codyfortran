! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> FFTW3-backed implementation of the SoftYukawa interaction on a linear grid.
!>
!> Implements convolution with the SoftYukawa kernel using FFTs for improved
!> performance on larger grids. Requires FFTW3 to be available and initialized
!> by the build/system configuration.
module M_SysInteraction_Linear_SoftYukawa_Fftw
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the FFTW3-based SoftYukawa implementation to the core interaction API.
    module subroutine SysInteraction_Linear_SoftYukawa_Fftw_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
