! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Linear.f90
!> @brief Linear-grid kinetic operators T̂ = −(1/2m)∇² on Cartesian grids.
!>
!> @details
!> This module represents the **linear (Cartesian, uniform 1D)** grid branch of
!> the kinetic subsystem. It provides the Laplacian operator ∇² on uniform grids
!> using either finite-difference stencils or FFT-based spectral methods.
!>
!> Physics
!> -------
!> For a particle of mass m on a uniform grid, the kinetic energy operator is:
!>
!>     T̂ψ(x) = −(1/2m) d²ψ/dx²
!>
!> This is discretized either by:
!>   - **FinDiff**: Central finite-difference stencil (O(h²) or higher)
!>   - **Fourier**: FFT → multiply by −k² → IFFT (spectral accuracy)
!>
!> The Fourier method assumes periodic boundary conditions; the FinDiff method
!> uses zero (Dirichlet) boundary conditions at the grid ends.
!>
!> @see M_SysKinetic_Linear_Laplacian for backend selection
module M_SysKinetic_Linear
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Fabricate the linear-grid kinetic backend.
    !>
    !> Reads configuration and assigns the setup pointer and operator routines.
    !> Currently supports: `laplacian` (FinDiff or Fourier).
    module subroutine SysKinetic_Linear_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to setup for the linear-grid kinetic backend.
  !>
  !> Precomputes FD stencils or spectral k-vectors, allocates workspace, and
  !> validates alignment/contiguity assumptions for the selected method.
  procedure(I_SysKinetic_Linear_Setup), pointer :: SysKinetic_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize the linear-grid kinetic backend (allocate/precompute data).
    subroutine I_SysKinetic_Linear_Setup
    end subroutine
  end interface

end module
