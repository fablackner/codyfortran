! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief FEM-grid (grid.fem) interaction interface.
!>
!> @details Pair densities and interaction potentials are nodal arrays on the
!> interior nodes of the spectral-element mesh (Grid_nPoints entries). The
!> family binds the parts shared by all nodal kernels:
!>   - FillInteractionSrc: the weighted source src = conjg(a) b w, which is
!>     exactly the Galerkin load of the pair density rho = conjg(a) b
!>     (Grid_Fem_SolvePoisson takes it as is);
!>   - MultiplyWithInteractionPotential: point-wise V orb;
!>   - ConjugateInteractionPotential: conjg(V) (nodal storage).
!> The model binds FillInteractionPotential.
!>
!> **Available implementations:**
!>   - `Coulomb`: 1/|r1 - r2| via the Galerkin Poisson solve with free-space
!>     (multipole) box values
module M_SysInteraction_Fem
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Check grid.fem, bind the nodal source/multiply/conjugate procedures and
    !> dispatch to the interaction model.
    module subroutine SysInteraction_Fem_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
