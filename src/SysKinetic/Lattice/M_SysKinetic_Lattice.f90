! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Lattice-based kinetic operators (tight-binding style) and their setup.
!>
!> This module hosts the lattice branch of the kinetic subsystem. It exposes a
!> `SysKinetic_Setup` pointer (assigned by the lattice fabric) and is intended
!> to bind concrete implementations such as nearest-neighbor hopping on linear
!> lattices or more elaborate stencils.
!>
!> Usage
!> 1) Call `SysKinetic_Lattice_Fabricate` to select the concrete lattice scheme
!>    and assign `SysKinetic_Setup`.
!> 2) Call `SysKinetic_Setup` prior to time propagation to allocate coefficients
!>    and validate dimensions.
module M_SysKinetic_Lattice
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the lattice kinetic backend and assign procedure pointers.
    !>
    !> Reads configuration (e.g., hopping amplitudes, connectivity, boundary
    !> conditions) and binds the setup pointer to the chosen scheme (e.g.,
    !> nearest-neighbor on a linear lattice). Must run before `SysKinetic_Setup`.
    module subroutine SysKinetic_Lattice_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup procedure for initializing the lattice kinetic system.
  !>
  !> Precomputes lattice-specific data such as neighbor lists, hopping matrices
  !> or direction-dependent coefficients required during time stepping.
  procedure(I_SysKinetic_Lattice_Setup), pointer :: SysKinetic_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize the lattice kinetic backend (allocate and precompute data).
    subroutine I_SysKinetic_Lattice_Setup
    end subroutine
  end interface

end module
