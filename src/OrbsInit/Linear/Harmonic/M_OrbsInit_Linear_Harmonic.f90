! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Quantum harmonic oscillator eigenstate initializer for 1D linear grids.
!>
!> @details
!> Implements orbital initialization using the exact eigenstates of the 1D
!> quantum harmonic oscillator:
!>
!>   ψ_n(x) ∝ H_n(ξ) · exp(-ξ²/2)
!>
!> where:
!> - ξ = (x - x₀) / a  is the dimensionless coordinate
!> - a = 1/√ω         is the oscillator length scale
!> - H_n(ξ)           is the physicist's Hermite polynomial of degree n
!> - x₀               is the trap center position
!> - ω                is the angular frequency
!>
!> The orbital index maps to quantum number as: n = index - 1
!> (so index=1 gives the ground state n=0, index=2 gives n=1, etc.)
!>
!> JSON Configuration
!> ------------------
!>   "orbsInit": {
!>     "linear": {
!>       "harmonic": {
!>         "position": 0.0,    // trap center x₀ (default: 0.0)
!>         "omega": 1.0        // angular frequency ω (default: 1.0)
!>       }
!>     }
!>   }
!>
!> @note Normalization is performed by the parent Linear module after sampling.
module M_OrbsInit_Linear_Harmonic
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate a harmonic-oscillator based linear initializer.
    !> Configures parameters for harmonic orbitals on a 1D grid and wires the
    !> callouts in the linear backend (e.g., ground/excited states centered at a
    !> specified position with frequency `omega`).
    module subroutine OrbsInit_Linear_Harmonic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Center position x₀ of the harmonic trap (same units as grid coordinate).
  real(R64) :: OrbsInit_Linear_Harmonic_position

  !> Angular frequency ω of the harmonic oscillator.
  !> Determines the oscillator length scale a = 1/√ω.
  real(R64) :: OrbsInit_Linear_Harmonic_omega

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
