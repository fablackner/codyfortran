! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Central interaction API for the code base.
!>
!> This module exposes a small, uniform set of procedure pointers that higher
!> layers call to construct and apply interaction potentials. At runtime,
!> `SysInteraction_Fabricate` binds these pointers to concrete implementations
!> (e.g., lattice, linear, or ylm back-ends) based on the chosen configuration.
!>
!> No interaction logic lives here directly; instead this module defines the
!> abstract interfaces and the public pointers that are wired up during
!> fabrication.
module M_SysInteraction
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Runtime factory that reads the user configuration and binds all
    !> `SysInteraction_*` procedure pointers to a concrete implementation.
    !>
    !> Typical responsibilities:
    !> - parse input parameters (e.g. from JSON)
    !> - select the grid/domain back-end (lattice/linear/ylm)
    !> - select the specific interaction model (e.g. on-site, Yukawa, Coulomb)
    !> - set feature flags like `SysInteraction_timeIndependentQ`
    module subroutine SysInteraction_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> True if the interaction is independent of body types (species). When set,
  !> callers may cache results across body-type combinations.
  logical :: SysInteraction_bodyTypeIndependentQ = .false.

  !> True if the interaction is time independent. When set, time arguments in
  !> back-end routines can be ignored and results cached/reused.
  logical :: SysInteraction_timeIndependentQ = .false.

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to a cheap, idempotent setup routine for the selected interaction
  !> back-end. Use this to allocate caches, precompute kernels, etc.
  procedure(I_SysInteraction_Setup), pointer :: SysInteraction_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initializes the selected interaction implementation (allocate tables,
    !> FFT plans, masks, boundary helpers, ...). Should be safe to call
    !> multiple times.
    subroutine I_SysInteraction_Setup
    end subroutine
  end interface

  !> Pointer to a routine that builds the full interaction potential resulting
  !> from the provided source term (accumulated over all source grid points).
  procedure(I_SysInteraction_FillInteractionPotential), pointer :: SysInteraction_FillInteractionPotential

  abstract interface
    subroutine I_SysInteraction_FillInteractionPotential(interactionPotential, src, time, bt1_, bt2_)
      import :: I32, R64
      !> Output: accumulated interaction potential on the active grid.
      complex(R64), intent(out), allocatable :: interactionPotential(:)
      !> Input: complex source term sampled on grid points.
      complex(R64), intent(in), contiguous  :: src(:)
      !> Physical time associated with the source/potential.
      real(R64), intent(in)     :: time
      !> Optional body type (species) index for the potential target.
      integer(I32), intent(in), optional :: bt1_
      !> Optional body type (species) index for the source.
      integer(I32), intent(in), optional :: bt2_
    end subroutine
  end interface

  !> Pointer to a routine that constructs a source term from a pair of orbitals
  !> (typically a bra-ket product) on the active grid.
  procedure(I_SysInteraction_FillInteractionSrc), pointer :: SysInteraction_FillInteractionSrc
  abstract interface
    !> Build a complex source term from two orbitals (e.g., rho = conjg(psi_i)*psi_j).
    subroutine I_SysInteraction_FillInteractionSrc(src, orbConjg, orb)
      import :: R64
      !> Output: complex source term over grid points.
      complex(R64), intent(out), allocatable :: src(:)
      !> Input: conjugated orbital (bra).
      complex(R64), intent(in), contiguous :: orbConjg(:)
      !> Input: orbital (ket).
      complex(R64), intent(in), contiguous :: orb(:)
    end subroutine
  end interface

  !> Pointer to a routine that applies the interaction potential to an orbital
  !> (point-wise multiplication on the active grid).
  procedure(I_SysInteraction_MultiplyWithInteractionPotential), pointer :: SysInteraction_MultiplyWithInteractionPotential
  abstract interface
    !> Multiply `orb` by `interactionPotential` element-wise and return the result.
    subroutine I_SysInteraction_MultiplyWithInteractionPotential(resOrb, interactionPotential, orb)
      import :: R64
      !> Output: resulting orbital after multiplication.
      complex(R64), intent(out), contiguous :: resOrb(:)
      !> Input: interaction potential over grid points.
      complex(R64), intent(in), contiguous  :: interactionPotential(:)
      !> Input: orbital to be multiplied.
      complex(R64), intent(in), contiguous  :: orb(:)
    end subroutine
  end interface

end module
