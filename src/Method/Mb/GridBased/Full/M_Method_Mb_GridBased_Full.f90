! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Full-expansion, grid-based representation of many-body wavefunctions.
!>
!> This variant represents the many-body state directly on the spatial grid (no
!> orbital/determinant expansions). It targets exact-time propagation within the
!> chosen discretization and enforces proper (anti-)symmetry based on statistics.
!>
!> State layout
!> - The packed state typically stores the full N-body wavefunction sampled on
!>   the direct-product grid. The ordering of grid indices follows the global
!>   body-ordering (compatible with `Method_Mb_*` index ranges).
!>
!> Operators
!> - Kinetic, potential, and interaction applications are provided by the
!>   grid-based operator interfaces and bound by the fabricate routine. They
!>   operate on 1D/2D slices extracted from the full state with statistics-aware
!>   symmetrization/antisymmetrization handled by the method.
module M_Method_Mb_GridBased_Full
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds Full-specific procedure pointers (setup, energy, time
    !> derivative, and grid-operator hooks) at runtime and initializes the
    !> direct grid representation according to the configuration (grid sizes,
    !> boundary conditions, statistics).
    module subroutine Method_Mb_GridBased_Full_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
