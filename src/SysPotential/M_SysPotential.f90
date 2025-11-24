! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Grid-agnostic external potential interface.
!>
!> This module defines the abstract, runtime-wired interface that all external
!> potentials in CodyFortran conform to. It holds feature flags and procedure
!> pointers which are assigned at runtime by `SysPotential_Fabricate`, based on
!> the selected backend (e.g., Linear, Lattice, Ylm) and the JSON configuration.
!>
!> No concrete physics is implemented here; only the contract that concrete
!> backends must satisfy and the pointers through which the rest of the code
!> interacts with an external potential.
module M_SysPotential
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Build-time independent factory that wires the potential at runtime.
    !>
    !> The factory inspects the runtime configuration (JSON) and assigns the
    !> procedure pointers exported by this module to concrete implementations
    !> provided by the chosen backend. It may also parse global parameters and
    !> set the optimization flags below.
    module subroutine SysPotential_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> True if the potential does not depend on the particle/body type.
  !> When set, callers may compute and reuse a single potential for all types.
  logical :: SysPotential_bodyTypeIndependentQ = .false.

  !> True if the potential is time-independent.
  !> When set, callers may cache the external potential across time steps.
  logical :: SysPotential_timeIndependentQ = .false.

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Setup hook for the potential backend (allocated buffers, masks, caches).
  !>
  !> Called once after fabrication and whenever a reinitialization is required
  !> due to configuration changes. No grid data is passed here; backends fetch
  !> any required global state from their owning grid modules.
  procedure(I_SysPotential_Setup), pointer :: SysPotential_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize the potential backend.
    subroutine I_SysPotential_Setup
    end subroutine
  end interface

  !> Compute V(r,t) on the active grid.
  !>
  !> The callee allocates and fills `externalPotential` with the complex values
  !> of the external potential at the active grid resolution. If the potential
  !> is diagonal in the chosen basis this is simply the point-wise potential
  !> V(r,t); more general operators may still expose an effective diagonal form
  !> expected by the caller.
  procedure(I_SysPotential_FillExternalPotential), pointer :: SysPotential_FillExternalPotential
  abstract interface
    !> Fill the external potential array for a given time (and optional body type).
    subroutine I_SysPotential_FillExternalPotential(externalPotential, time, bt_)
      import :: I32, R64
      !> Output: allocated array holding V on the active grid (shape defined by backend).
      complex(R64), intent(out), allocatable :: externalPotential(:)
      !> Current simulation time for which the potential is evaluated.
      real(R64), intent(in)  :: time
      !> Optional body type selector; ignored when SysPotential_bodyTypeIndependentQ is true.
      integer(I32), intent(in), optional :: bt_
    end subroutine
  end interface

  !> Apply the diagonal external potential to an orbital: dOrb := V * orb.
  !>
  !> This is a point-wise multiplication (Hadamard product) in the representation
  !> where the potential is diagonal. The arrays must be contiguous and conforming.
  procedure(I_SysPotential_MultiplyWithExternalPotential), pointer :: SysPotential_MultiplyWithExternalPotential
  abstract interface
    !> Multiply an input orbital with the external potential array.
    subroutine I_SysPotential_MultiplyWithExternalPotential(dOrb, externalPotential, orb)
      import :: R64
      !> Output: result orbital (same shape as inputs).
      complex(R64), intent(out), contiguous :: dOrb(:)
      !> Input: external potential array V.
      complex(R64), intent(in), contiguous :: externalPotential(:)
      !> Input: orbital to be multiplied by V.
      complex(R64), intent(in), contiguous :: orb(:)
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
