! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Central facade for kinetic-energy operators used throughout the codebase.
!>
!> This module declares procedure pointers that are bound at runtime by
!> `SysKinetic_Fabricate`, based on the simulation setup (grid type, numerical
!> scheme, masses, etc.). Callers use these pointers without depending on a
!> specific backend (linear grid, lattice, ylm, finite-difference, Fourier,…).
!>
!> Contract
!> - Fabrication selects and assigns concrete procedures before first use.
!> - Setup is idempotent and prepares any masks, work arrays, or constants.
!> - Kinetic application operates on contiguous complex arrays representing a
!>   single orbital; an optional body-type index allows mass scaling.
module M_SysKinetic
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Build-time independent factory for kinetic operators.
    !>
    !> Reads user configuration (e.g., JSON input) and assigns the procedure
    !> pointers exported by this module to concrete implementations, such as
    !> linear-grid Laplacian (finite-difference/Fourier), lattice hopping, or
    !> spherical-ylm radial operators. Must be invoked exactly once per
    !> simulation context before calling any other procedure in this module.
    module subroutine SysKinetic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> True if the kinetic operator does not depend on time
  logical :: SysKinetic_timeIndependentQ = .false.

  !> True if the kinetic operator is body-type independent
  logical :: SysKinetic_bodyTypeIndependentQ = .false.

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup procedure for initializing the kinetic subsystem.
  !>
  !> Responsibilities
  !> - Allocate/capture grid-dependent constants and masks
  !> - Precompute operator coefficients (e.g., k^2 arrays for Fourier, FD stencils)
  !> - Validate dimensional compatibility with the chosen backend
  procedure(I_SysKinetic_Setup), pointer :: SysKinetic_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initializes the kinetic subsystem, including any mask functions and
    !> numerical parameters needed by the selected kinetic implementation.
    subroutine I_SysKinetic_Setup
    end subroutine
  end interface

  !> Pointer to the procedure that applies the kinetic operator to a single orbital.
  !>
  !> Notes
  !> - The concrete operator depends on the selected backend (e.g., −(1/2m)∇² on
  !>   a linear grid, tight-binding hopping on a lattice, or the radial part in ylm).
  !> - Arrays must be contiguous; shapes and boundary conditions are backend-specific.
  !> - If provided, `bt_` selects a body-type specific mass for the 1/(2 m_bt) factor.
  procedure(I_SysKinetic_MultiplyWithKineticOp), pointer :: SysKinetic_MultiplyWithKineticOp
  abstract interface
    !> Applies the kinetic-energy operator to the input orbital and returns the result.
    !>
    !> Parameters
    !> - dOrb: output orbital after applying the operator (contiguous)
    !> - orb: input orbital to be transformed (contiguous)
    !> - time: simulation time for time-dependent operators (if any)
    !> - bt_: optional body-type index for mass scaling
    subroutine I_SysKinetic_MultiplyWithKineticOp(dOrb, orb, time, bt_)
      import :: I32, R64
      !> Output orbital after applying the kinetic operator.
      complex(R64), intent(out), contiguous :: dOrb(:)
      !> Input orbital.
      complex(R64), intent(in), contiguous :: orb(:)
      !> Current time, allowing for time-dependence.
      real(R64), intent(in) :: time
      !> Body type index for mass-scaling.
      integer(I32), intent(in), optional :: bt_
    end subroutine
  end interface

end module
