! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_TwoRdmInit_Excited.f90
!> @brief Excited-state specialization of the 2-RDM initializer.
!>
!> @details
!> This module provides a 2-RDM initializer that constructs the two-body
!> reduced density matrix from products of 1-RDMs. This is appropriate for:
!>
!>   - **Hartree-Fock initial states**: Single Slater determinants where
!>     D^2(i,j;k,l) = D^1(i,k) D^1(j,l) - D^1(i,l) D^1(j,k)
!>
!>   - **Single excitations**: Applying creation/annihilation operators
!>     to modify orbital occupations from a reference configuration.
!>
!> @par Physical Interpretation
!> For a state |Phi> = prod_i a^dag_i |0>, the 2-RDM factorizes:
!>
!>   Same-spin:     D^2_aa(ij;kl) = n_i n_j (delta_ik delta_jl - delta_il delta_jk)
!>   Opposite-spin: D^2_ab(ij;kl) = n_i n_j delta_ik delta_jl
!>
!> where n_i are orbital occupation numbers (0 or 1 for fermions).
!>
!> @par JSON Configuration
!> @code{.json}
!> {
!>   "twoRdmInit": {
!>     "excited": {
!>       "creates": [3],        // Orbitals to populate (optional, default: [1])
!>       "destroys": [1],       // Orbitals to depopulate (optional, default: [1])
!>       "bodyType1": 1,        // Body type for creates/destroys (optional)
!>       "bodyType2": 2         // Body type for second species (optional)
!>     }
!>   }
!> }
!> @endcode
!>
!> @note The creates/destroys parameters are read but currently not applied
!>       in a fully general way. The implementation builds diagonal 1-RDMs.
module M_TwoRdmInit_Excited
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Bind `M_TwoRdmInit` procedure pointers to the excited-state
    !>        implementation and read configuration parameters.
    !>
    !> @details
    !> Reads from `"twoRdmInit.excited"`:
    !>   - `creates`: Array of orbital indices to occupy (default: [1])
    !>   - `destroys`: Array of orbital indices to vacate (default: [1])
    !>   - `bodyType1`, `bodyType2`: Spin species indices (default: 1, 2)
    !>
    !> Binds `TwoRdmInit_InitializeBlock` -> `InitializeBlock`
    module subroutine TwoRdmInit_Excited_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

end module
