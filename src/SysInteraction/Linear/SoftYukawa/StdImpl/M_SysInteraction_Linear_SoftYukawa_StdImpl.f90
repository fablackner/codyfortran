! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Direct convolution implementation of SoftYukawa interaction.
!>
!> @details Uses the `ConvolutionIntegral` utility for O(N²) real-space
!> convolution. Exact and dependency-free, but slow for large grids (N > 1000).
!> Prefer the FFT variant for production runs on large systems.
module M_SysInteraction_Linear_SoftYukawa_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the default linear SoftYukawa implementation to the core interaction
    !> API exposed by `M_SysInteraction`.
    module subroutine SysInteraction_Linear_SoftYukawa_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
