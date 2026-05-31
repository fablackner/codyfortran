! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Reference (Green's function) implementation of radial Coulomb solver.
!>
!> @details Computes the radial potential using the analytical Green's function:
!>    Vₗₘ(r) = (4π/(2l+1)) ∫ (r<ˡ/r>ˡ⁺¹) ρₗₘ(r') dr'
!>
!> **Complexity:** O(N²) where N is the number of radial points.
!> **Use case:** Reference/validation, small grids (N < 100).
module M_SysInteraction_Ylm_Coulomb_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the default Ylm Coulomb solver to the core interaction API.
    module subroutine SysInteraction_Ylm_Coulomb_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
