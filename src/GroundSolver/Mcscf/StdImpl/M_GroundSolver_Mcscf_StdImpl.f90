! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Standard MCSCF implementation interface module.
!>
!> @details
!> Provides the factory routine that binds the MCSCF callbacks to a reference
!> implementation suitable for general grids. This backend:
!>
!>   - Operates on the packed MCTDHX-style state (CI coefficients ⊕ orbitals)
!>   - Uses DiagonalizerList(1) over the full CI space (dim = Coeffs_nCoeffs)
!>   - Uses DiagonalizerList(2) over the orbital space (dim = Grid_nPoints)
!>     with the RDM-weighted Fock operator
!>
!> ## JSON Configuration
!>
!> ```json
!> {
!>     "groundSolver": {
!>         "mcscf": {
!>             "stdImpl": { }
!>         }
!>     }
!> }
!> ```
!>
!> The `alpha` argument of `GroundSolver_Approach` mixes both the CI
!> coefficients with the CI ground eigenvector and the orbitals with the
!> Fock eigenvectors.
!>
!> @see M_GroundSolver_Mcscf, M_GroundSolver_Scf_StdImpl
module M_GroundSolver_Mcscf_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Bind MCSCF procedure pointers to the standard implementation.
    !>
    !> @details
    !> Assigns the following procedure pointers:
    !>   - `GroundSolver_Setup` → local `Setup`
    !>   - `GroundSolver_Approach` → local `Approach`
    !>   - `GroundSolver_Mcscf_HamiltonianAction` → local `HamiltonianAction`
    !>   - `GroundSolver_Mcscf_FockAction` → local `FockAction`
    !>
    !> @pre JSON key `groundSolver.mcscf.stdImpl` exists
    module subroutine GroundSolver_Mcscf_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
