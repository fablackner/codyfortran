! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default real-space implementation of the SoftYukawa interaction on a linear grid.
!>
!> Provides a simple, dependency-free realization (e.g., direct convolution) of
!> the SoftYukawa kernel. For large grids, FFT-backed variants may be preferred.
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
