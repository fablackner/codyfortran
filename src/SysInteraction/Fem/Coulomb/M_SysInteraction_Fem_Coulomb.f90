! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Coulomb interaction strength/|r1 - r2| on grid.fem.
!>
!> @details The pair potential of the weighted source src = conjg(a) b w is
!>   V = strength * Grid_Fem_SolvePoisson(src),
!> the Galerkin solution of -lap V = 4 pi rho, rho = conjg(a) b, with
!> free-space box values from the multipole expansion of the source
!> (grid.fem.multipoleLmax). Real and imaginary parts are solved separately
!> (a real source costs one FDM solve, a complex one two). The kernel and the
!> weights are real, so the potential of the swapped pair is the conjugate
!> (SysInteraction_conjSymmetricQ); time and body types do not enter
!> (Method's bt1_ = bt2_ = -1 are ignored).
!>
!> JSON Configuration
!> ------------------
!>   "sysInteraction": { "fem": { "coulomb": { "strength": 1.0 } } }
module M_SysInteraction_Fem_Coulomb
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Read the strength, set the flags and bind FillInteractionPotential.
    module subroutine SysInteraction_Fem_Coulomb_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Coupling strength multiplying 1/|r1 - r2| (1.0 for electrons).
  real(R64) :: SysInteraction_Fem_Coulomb_strength = 1.0_R64

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
