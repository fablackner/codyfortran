! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Ylm-optimized TDHx implementation interface module.
!>
!> @details
!> Provides the factory routine and interface for a spherical-harmonics (Ylm)
!> optimized TDHx ground-state solver. This variant exploits angular momentum
!> symmetry to reduce computational cost:
!>
!> **Key Optimization**: Instead of diagonalizing a single Grid_nPoints×Grid_nPoints
!> Fock matrix, this backend diagonalizes (lmax+1) smaller nRadial×nRadial matrices,
!> one per angular momentum channel l. This is valid for spherically symmetric
!> Hamiltonians (e.g., atoms) where orbitals factorize as φ_{nlm}(r) = R_nl(r)·Y_lm(θ,φ).
!>
!> ## Usage
!>
!> Requires a DiagonalizerList with (lmax+1) entries, one per l-channel. The user
!> must provide wrapper callbacks that fix the l-value for each diagonalizer:
!>
!> ```fortran
!> DiagonalizerListInput(1) % ApplyMatOnVec => ApplyMatOnVecL0  ! l=0
!> DiagonalizerListInput(2) % ApplyMatOnVec => ApplyMatOnVecL1  ! l=1
!> ! ...
!> ```
!>
!> @see M_GroundSolver_Tdhx_StdImpl for the general-grid alternative
module M_GroundSolver_Tdhx_YlmOpt
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Bind Ylm-optimized TDHx callbacks for the ground-state solver.
    !>
    !> @details
    !> Assigns the following procedure pointers:
    !>   - `GroundSolver_Setup` → local `Setup`
    !>   - `GroundSolver_Approach` → local `Approach`
    !>   - `GroundSolver_Tdhx_YlmOpt_HartreeFockAction` → local `HartreeFockAction`
    !>
    !> @pre JSON key `groundSolver.tdhx.ylmOpt` exists
    !> @pre Grid is configured as Ylm (spherical harmonics)
    module subroutine GroundSolver_Tdhx_YlmOpt_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to the Ylm radial Fock operator action for a single l-channel.
  !>
  !> Used as the matvec callback for per-l diagonalizers. The caller must wrap
  !> this in l-specific callbacks (one per diagonalizer).
  procedure(I_GroundSolver_Tdhx_YlmOpt_HartreeFockRadialAction), pointer :: GroundSolver_Tdhx_YlmOpt_HartreeFockAction
  abstract interface
    !> @brief Apply the radial Fock operator for angular momentum channel l.
    !>
    !> @details
    !> Computes dOrbLm = F̂_l · orbLm where F̂_l is the radial Fock operator
    !> projected onto channel l. This includes:
    !>   - Radial kinetic: T̂_l = -½d²/dr² + l(l+1)/(2r²)
    !>   - Radial external potential: V̂_ext(r)
    !>   - Radial mean-field potential: precomputed monopole (l=0) projection
    !>   - Exchange: computed on-the-fly with full angular coupling
    !>
    !> @param[out] dOrbLm  Result of F̂_l applied to orbLm, dimension(nRadial)
    !> @param[in]  orbLm   Input radial orbital for channel l
    !> @param[in]  l       Angular momentum quantum number
    !> @param[in]  time    Time parameter (typically 0 for ground state)
    subroutine I_GroundSolver_Tdhx_YlmOpt_HartreeFockRadialAction(dOrbLm, orbLm, l, time)
      import :: I32, R64
      !> Output: radial Fock operator applied to the orbital, F̂_l·orbLm
      complex(R64), intent(out), contiguous, target :: dOrbLm(:)
      !> Input radial orbital for channel l
      complex(R64), intent(in), contiguous, target :: orbLm(:)
      !> Angular momentum quantum number (l = 0, 1, ..., lmax)
      integer(I32), intent(in) :: l
      !> Time parameter (passed through to potential routines)
      real(R64), intent(in) :: time
    end subroutine
  end interface

end module

