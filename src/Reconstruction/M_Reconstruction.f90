! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Core reconstruction wiring for RDM-based dynamics.
!>
!> @details
!> **Purpose**: This module provides the public interface for BBGKY hierarchy
!> closure in time-dependent 2-RDM (TD-2RDM) simulations. The equation of motion
!> for the 2-RDM involves the 3-RDM through the collision operator:
!>
!>     d(D²)/dt = ... + Tr₃[Ŵ, D³]
!>
!> Reconstruction schemes approximate D³ from lower-order RDMs (D¹, D²), avoiding
!> the exponential cost of tracking the full 3-RDM.
!>
!> **Design**: This module does not implement any reconstruction scheme itself.
!> Instead, it declares procedure pointers that the rest of the code calls, and
!> provides a factory routine to bind concrete implementations at runtime. The
!> JSON configuration selects among available families:
!>
!> - `generic`   : Spin-unresolved, full tensor representation
!> - `singlet`   : Spin-singlet (S=0), exploits α-β block structure
!> - `singletTr0`: Trace-free singlet representation with explicit 1-RDM
!>
!> **Usage pattern**:
!> 1. Call `Reconstruction_Fabricate` once during program setup to bind the
!>    procedure pointers based on JSON configuration.
!> 2. Call `Reconstruction_Setup` (if the chosen scheme requires it) to
!>    allocate workspace and read scheme-specific options.
!> 3. During time propagation, call `Reconstruction_AddCollisionOp` to
!>    accumulate the collision-operator contribution to dD²/dt.
!>
!> @see M_Reconstruction_Generic, M_Reconstruction_Singlet, M_Reconstruction_SingletTr0
module M_Reconstruction
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Attach concrete reconstruction procedures for the selected scheme.
    !>
    !> Called once during initialization. The implementation chooses a
    !> reconstruction family (generic, singlet, etc.), wires the procedure
    !> pointers exported by this module to the corresponding concrete routines,
    !> and may perform lightweight global initialization. Any heavy
    !> allocation/initialization should be deferred to `Reconstruction_Setup`.
    module subroutine Reconstruction_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to the optional setup routine for the chosen scheme.
  !>
  !> Defaults to a no-op. Schemes requiring workspace allocation or option
  !> parsing override this during fabrication.
  procedure(I_Reconstruction_Setup), pointer :: Reconstruction_Setup => NoOpProcedures_Setup
  abstract interface
    !> @brief Initialize the chosen reconstruction scheme.
    !>
    !> @details
    !> Implementations typically perform one or more of:
    !> - Allocate persistent work arrays (e.g., rdm3, collOp2)
    !> - Precompute index maps or symmetry factors
    !> - Read scheme-specific JSON options
    !> - Load ML model files for neural-network schemes
    !>
    !> Schemes that do not require setup may leave the default no-op.
    subroutine I_Reconstruction_Setup
    end subroutine
  end interface

  !> @brief Pointer to the routine that accumulates the 2-body collision operator.
  !>
  !> This is the main entry point called during each time derivative evaluation.
  procedure(I_Reconstruction_AddCollisionOp), pointer :: Reconstruction_AddCollisionOp
  abstract interface
    !> @brief Add the collision-operator contribution to d(D²)/dt in-place.
    !>
    !> @details
    !> Computes the BBGKY collision term Tr₃[Ŵ, D³] where D³ is reconstructed
    !> from D¹ and D² using the active scheme (Valdemoro, NY, neural network, etc.).
    !>
    !> **Important**: This routine ADDS to `dTwoRdmState`; it does not overwrite.
    !> The tensor layouts follow project-wide conventions:
    !> - Generic: full nOrbs⁴ tensors
    !> - Singlet/SingletTr0: α-β block rdm2ab(nα,nα,nβ,nβ)
    !>
    !> @param[inout] dTwoRdmState Packed derivative vector; collision term added in-place
    !> @param[in]    twoRdmState  Packed 2-RDM state at current time
    !> @param[in]    h2           Two-body Hamiltonian h(p,q,r,s) = ⟨pq|Ŵ|rs⟩
    subroutine I_Reconstruction_AddCollisionOp(dTwoRdmState, twoRdmState, h2)
      import :: R64
      complex(R64), intent(inout), contiguous  :: dTwoRdmState(:)
      complex(R64), intent(in), contiguous     :: twoRdmState(:)
      complex(R64), intent(in), contiguous     :: h2(:, :, :, :)
    end subroutine
  end interface

end module

