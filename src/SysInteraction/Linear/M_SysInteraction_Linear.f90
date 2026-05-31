! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Linear-grid (1D real-space) interaction interface.
!>
!> @details This module registers interaction models for uniformly spaced 1D
!> grids where the interaction potential is computed via real-space convolution:
!>    V(x) = ∫ w(|x - x'|) ρ(x') dx'
!>
!> **Available implementations:**
!>   - `SoftYukawa`: Screened Coulomb/Yukawa kernel with softening
!>
!> The convolution can be computed via direct O(N²) integration (`StdImpl`)
!> or FFT-accelerated O(N log N) method (`Fftw3`).
module M_SysInteraction_Linear
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind linear-grid interaction procedures according to configuration.
    module subroutine SysInteraction_Linear_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
