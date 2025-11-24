! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Linear-grid kinetic operators and setup glue.
!>
!> This module represents the linear (Cartesian, uniform 1D/ND) grid branch of
!> the kinetic subsystem. It exposes a `SysKinetic_Setup` pointer and a
!> `SysKinetic_Linear_Fabricate` entry point that chooses a concrete Laplacian
!> backend (e.g., finite-difference or Fourier/spectral) and binds the pointers.
module M_SysKinetic_Linear
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the linear-grid kinetic backend.
    !>
    !> Reads configuration (grid spacing, boundary conditions, mass model, and
    !> method selection like FD/Fourier) and assigns the setup pointer and the
    !> operator application routines in the parent facade.
    module subroutine SysKinetic_Linear_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup procedure for initializing the linear-grid backend.
  !>
  !> Precomputes FD stencils or spectral k-vectors, allocates workspace, and
  !> validates alignment/contiguity assumptions for the selected method.
  procedure(I_SysKinetic_Linear_Setup), pointer :: SysKinetic_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize the linear-grid kinetic backend (allocate/precompute data).
    subroutine I_SysKinetic_Linear_Setup
    end subroutine
  end interface

end module
