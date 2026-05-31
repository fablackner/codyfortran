! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief FFT-accelerated implementation of SoftYukawa interaction.
!>
!> @details Uses FFTW3 for O(N log N) convolution via the convolution theorem:
!>    V = IFFT( FFT(kernel) * FFT(ρ) )
!>
!> Requires FFTW3 library. Recommended for grids with N > 500 points where
!> the O(N log N) scaling significantly outperforms direct O(N²) summation.
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
