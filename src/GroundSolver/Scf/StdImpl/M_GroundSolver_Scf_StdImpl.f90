! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Standard SCF implementation interface module.
!>
!> @details
!> Provides the factory routine that binds SCF callbacks to a reference
!> implementation that works on any grid back-end. This backend:
!>
!>   - Operates on the full Grid_nPoints-dimensional Hilbert space
!>   - Uses a single diagonalizer for all orbitals
!>   - Suitable for 1D (Linear), 2D and Ylm grids
!>
!> For spherical-harmonics (Ylm) grids with angular symmetry, the `ylmOpt`
!> backend is faster since it diagonalizes smaller per-l-channel matrices.
!>
!> @see M_GroundSolver_Scf_YlmOpt
module M_GroundSolver_Scf_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Bind SCF procedure pointers to the standard implementation.
    !>
    !> @details
    !> Assigns the following procedure pointers:
    !>   - `GroundSolver_Setup` → local `Setup`
    !>   - `GroundSolver_Approach` → local `Approach`
    !>   - `GroundSolver_Scf_HartreeFockAction` → local `HartreeFockAction`
    !>
    !> @pre JSON key `groundSolver.scf.stdImpl` exists
    module subroutine GroundSolver_Scf_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module

