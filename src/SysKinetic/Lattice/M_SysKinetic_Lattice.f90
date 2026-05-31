! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Lattice.f90
!> @brief Lattice-based kinetic operators (tight-binding style).
!>
!> @details
!> This module hosts the **lattice branch** of the kinetic subsystem, implementing
!> kinetic energy as nearest-neighbor hopping on a discrete 3D lattice.
!>
!> Physics Background
!> ------------------
!> In tight-binding models, the kinetic energy is represented by hopping terms:
!>
!>     T̂ = −Σ_{⟨i,j⟩} t_{ij} (ĉ†_i ĉ_j + h.c.)
!>
!> where t_{ij} is the hopping amplitude between adjacent sites i and j. This
!> discrete Hamiltonian replaces the continuum Laplacian and is the foundation
!> for Hubbard-type models.
!>
!> Usage
!> -----
!> 1. Call `SysKinetic_Lattice_Fabricate` to select the concrete lattice scheme
!>    and assign `SysKinetic_Setup`.
!> 2. Call `SysKinetic_Setup` prior to time propagation to allocate coefficients
!>    and validate dimensions.
!>
!> JSON Configuration
!> ------------------
!> @code{.json}
!> {
!>   "sysKinetic": {
!>     "lattice": {
!>       "nearestNeighbor": {
!>         "hoppX": 1.0,
!>         "hoppY": 1.0,
!>         "hoppZ": 1.0
!>       }
!>     }
!>   }
!> }
!> @endcode
!>
!> @see M_SysKinetic_Lattice_NearestNeighbor, M_Grid_Lattice
module M_SysKinetic_Lattice
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Fabricate the lattice kinetic backend and assign procedure pointers.
    !>
    !> Reads configuration (hopping amplitudes, connectivity, boundary
    !> conditions) and binds the setup pointer to the chosen scheme.
    !> Must run before `SysKinetic_Setup`.
    module subroutine SysKinetic_Lattice_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to setup for the lattice kinetic subsystem.
  !>
  !> Precomputes lattice-specific data such as neighbor lists, hopping matrices,
  !> or direction-dependent coefficients required during time stepping.
  procedure(I_SysKinetic_Lattice_Setup), pointer :: SysKinetic_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize lattice kinetic backend (allocate/precompute data).
    subroutine I_SysKinetic_Lattice_Setup
    end subroutine
  end interface

end module
