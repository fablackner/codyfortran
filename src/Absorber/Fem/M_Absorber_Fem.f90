! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief FEM-grid (grid.fem) absorber family and factory.
!>
!> @details
!> Absorbers for the 3D Cartesian spectral-element grid: a mask on the
!> interior nodes that damps the orbitals in boundary layers along the six
!> faces of the box (the box surface itself carries homogeneous Dirichlet
!> conditions, so without an absorber outgoing waves reflect there).
!>
!> @par Available Implementations
!> - **cosinus**: product over the axes of cos^(1/order) layer profiles (see
!>   `M_Absorber_Fem_Cosinus`)
!>
!> @see M_Absorber, M_Absorber_Fem_Cosinus
module M_Absorber_Fem
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Install a FEM-grid absorber implementation and read input.
    !>
    !> @details
    !> Checks grid.fem and dispatches to:
    !>
    !> - **absorber.fem.cosinus**: cosine-profile boundary layers
    !>
    !> @post `Absorber_Setup` and `Absorber_ApplyAbsorber` are bound.
    module subroutine Absorber_Fem_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
