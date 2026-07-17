! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Grid-agnostic gauge-coupling interface (velocity gauge).
!>
!> This module defines the abstract, runtime-wired interface for the one-body
!> light-matter coupling term of the Hamiltonian in velocity gauge and dipole
!> approximation,
!>
!>     H_A = A(t)·p / m + A(t)² / (2 m)
!>
!> which complements the multiplicative external potential (SysPotential) and
!> the kinetic operator (SysKinetic):
!>
!>     H = T + V_ext + H_A + W
!>
!> Unlike SysPotential, the gauge coupling is not diagonal in position space
!> (it contains the momentum operator), so it follows the SysKinetic operator
!> contract: a single application routine instead of a Fill/Multiply pair.
!>
!> No concrete physics is implemented here; procedure pointers are assigned at
!> runtime by `SysGauge_Fabricate` based on the selected backend (e.g., Ylm)
!> and the JSON configuration. Configure only when a laser coupling in
!> velocity gauge is wanted; length-gauge couplings (E(t)·z) remain ordinary
!> multiplicative potentials and belong to SysPotential.
module M_SysGauge
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Build-time independent factory that wires the gauge coupling at runtime.
    !>
    !> Inspects the runtime configuration (JSON) and assigns the procedure
    !> pointers exported by this module to concrete implementations provided
    !> by the chosen backend. It may also parse global parameters and set the
    !> optimization flags below.
    module subroutine SysGauge_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> True if the gauge coupling does not depend on the particle/body type.
  logical :: SysGauge_bodyTypeIndependentQ = .false.

  !> True if the gauge coupling is time-independent. A physical vector
  !> potential A(t) is time-dependent, so backends set this to .false.;
  !> the flag exists for symmetry with the other Hamiltonian-term modules.
  logical :: SysGauge_timeIndependentQ = .false.

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Setup hook for the gauge backend (precomputed coordinates, caches).
  !>
  !> Called once after fabrication and after `Grid_Setup`, since backends
  !> capture grid data (radial coordinates, derivative contexts) here.
  procedure(I_SysGauge_Setup), pointer :: SysGauge_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize the gauge backend.
    subroutine I_SysGauge_Setup
    end subroutine
  end interface

  !> Apply the gauge coupling to an orbital: dOrb = H_A(t) · orb.
  !>
  !> The operator contains the momentum operator and is therefore not a
  !> point-wise multiplication; shapes and boundary conditions are
  !> backend-specific, matching the active grid representation.
  procedure(I_SysGauge_MultiplyWithGaugeOp), pointer :: SysGauge_MultiplyWithGaugeOp
  abstract interface
    !> Apply the velocity-gauge coupling operator: dOrb = H_A · orb.
    subroutine I_SysGauge_MultiplyWithGaugeOp(dOrb, orb, time, bt_)
      import :: I32, R64
      !> Output orbital after applying the gauge coupling: dOrb = H_A · orb.
      complex(R64), intent(out), contiguous :: dOrb(:)
      !> Input orbital wavefunction on the spatial grid.
      complex(R64), intent(in), contiguous :: orb(:)
      !> Current simulation time at which A(t) is evaluated.
      real(R64), intent(in) :: time
      !> Optional body type index (1-based) selecting the mass in A·p/m.
      integer(I32), intent(in), optional :: bt_
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
