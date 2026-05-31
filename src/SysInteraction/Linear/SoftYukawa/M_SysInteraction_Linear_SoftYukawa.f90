! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Softened Yukawa/screened-Coulomb interaction on linear grids.
!>
!> @details Implements a flexible interaction kernel of the form:
!>    w(r) = strength * exp(-dampening * r) / (sqrt(r² + softening1²) + softening2)
!>
!> This kernel interpolates between several physically relevant limits:
!>   - Pure Coulomb: dampening=0, softening1→0, softening2=0 gives 1/r
!>   - Yukawa: dampening>0, softening→0 gives exp(-κr)/r (screened)
!>   - Soft-core: softening1>0 regularizes the r→0 singularity
!>   - Contact: large dampening localizes interaction to small r
!>
!> **JSON configuration:**
!> ```json
!> "sysInteraction": {
!>   "linear": {
!>     "softYukawa": {
!>       "strength": 1.0,
!>       "softening1": 1.0,
!>       "softening2": 0.0,
!>       "dampening": 0.0,
!>       "stdImpl": {}
!>     }
!>   }
!> }
!> ```
!>
!> **Implementations:**
!>   - `stdImpl`: Direct convolution O(N²), exact but slow for large grids
!>   - `fftw`: FFT-based convolution O(N log N), requires FFTW3
module M_SysInteraction_Linear_SoftYukawa
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Register the SoftYukawa model and expose its kernel to the linear-grid
    !> back-end (either direct or FFT-based).
    module subroutine SysInteraction_Linear_SoftYukawa_Fabricate
    end subroutine

    !> Kernel function K(r) of the SoftYukawa model evaluated at a scalar
    !> distance `r` (in the active grid's length units). The exact form is
    !> controlled by the parameters below.
    pure module function SysInteraction_Linear_SoftYukawa_Interaction(distance) result(res)
      real(R64) :: res
      real(R64), intent(in) :: distance
    end function
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Overall coupling strength of the kernel.
  real(R64) :: SysInteraction_Linear_SoftYukawa_strength
  !> Inner softening parameter (near-field regularization).
  real(R64) :: SysInteraction_Linear_SoftYukawa_softening1
  !> Outer softening/shape parameter.
  real(R64) :: SysInteraction_Linear_SoftYukawa_softening2
  !> Exponential damping coefficient controlling long-range decay.
  real(R64) :: SysInteraction_Linear_SoftYukawa_dampening

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
