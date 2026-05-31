! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Linear_Laplacian_Fourier.f90
!> @brief Fourier/spectral Laplacian for linear grids.
!>
!> @details
!> Implements a spectral ∇² using FFT-based methods on uniform grids. The
!> Laplacian is computed by:
!>
!>     ∇²ψ = FFT⁻¹[ −k² · FFT[ψ] ]
!>
!> where k is the wave vector corresponding to each Fourier mode.
!>
!> Boundary Conditions
!> -------------------
!> Assumes **periodic** boundary conditions: ψ(x_min) = ψ(x_max). For problems
!> with Dirichlet BCs, use the FinDiff backend instead.
!>
!> Advantages
!> ----------
!>   - Spectral accuracy (exponential convergence for smooth functions)
!>   - O(N log N) complexity via FFT
!>   - No dispersion errors (exact k² multiplication)
!>
!> Requirements
!> ------------
!> Requires `grid.linear.const` to be configured (uniform grid spacing).
!>
!> @see M_Utils_DerivativeFftw for the underlying FFT implementation
module M_SysKinetic_Linear_Laplacian_Fourier
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Configure the Fourier Laplacian and bind the operator.
    !>
    !> Validates that `grid.linear.const` is available, then assigns
    !> `SysKinetic_MultiplyWithKineticOp` and `SysKinetic_Setup`.
    module subroutine SysKinetic_Linear_Laplacian_Fourier_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
