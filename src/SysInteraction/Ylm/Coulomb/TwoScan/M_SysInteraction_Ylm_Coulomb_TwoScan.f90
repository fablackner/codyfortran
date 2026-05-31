! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Two-scan O(N) radial Coulomb solver.
!>
!> @details Exploits the separable structure of the Green's function to compute
!> the radial potential in two linear passes (prefix and suffix sums):
!>    C1(i) = Σⱼ≤ᵢ ρ(j) rⱼˡ        (forward scan)
!>    C2(i) = Σⱼ≥ᵢ ρ(j) / rⱼˡ⁺¹    (backward scan)
!>    V(i) = factor × (C1(i)/rᵢˡ⁺¹ + rᵢˡ × C2(i) - ρ(i)/rᵢ)
!>
!> **Complexity:** O(N) — optimal for large radial grids.
!> **Use case:** Production runs, large grids (N > 100).
module M_SysInteraction_Ylm_Coulomb_TwoScan
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the two-scan Coulomb solver to the core Ylm interaction hooks.
    module subroutine SysInteraction_Ylm_Coulomb_TwoScan_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
