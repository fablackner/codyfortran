! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> All-active fermionic configuration implementation.
!>
!> Concrete implementation of the all-active configuration element for fermions.
!> Provides the excitation application logic and statistics-specific setup and
!> fabrication steps. Fermionic operators obey anti-commutation relations and
!> enforce the Pauli exclusion principle; sign factors arise from permutation
!> parity when applying excitations.
module M_ConfigList_AllActive_Fermionic
  use M_Utils_Types
  use M_ConfigList, only: T_ConfigList_E
  use M_ConfigList_AllActive, only: T_ConfigList_E_AllActive

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocate a fermionic all-active configuration element from input at `path`.
    !>
    !> Returns the concrete element as a polymorphic
    !> `class(T_ConfigList_E)` to be stored in the global list.
    module subroutine ConfigList_E_AllActive_Fermionic_Allocate(e, path)
      class(T_ConfigList_E), allocatable, intent(out) :: e
      !> Path into the input document for this element's configuration.
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> Concrete all-active fermionic configuration element.
  !>
  !> Implements fermionic excitation logic (anti-commuting operators with sign)
  !> and provides the statistics-specific level-2 setup/fabrication hooks used
  !> by the base class.
  type, extends(T_ConfigList_E_AllActive) :: T_ConfigList_E_AllActive_Fermionic

  contains

    !> Apply creation and annihilation operators using fermionic anti-commutation rules.
    procedure :: ExciteConfiguration

    !> Fermionic level-2 setup invoked by the all-active base `Setup`.
    procedure :: SetupLevel2 => Setup

    !> Fermionic level-2 fabrication invoked by the all-active base `Fabricate`.
    procedure :: FabricateLevel2 => Fabricate

  end type

  interface
    !> Apply fermionic excitation operators to transform a configuration index.
    !>
    !> Fermionic operators obey the anti-commutation relations
    !> {a_i, a^\dagger_j} = \delta_{ij}, with {a_i, a_j} = {a^\dagger_i, a^\dagger_j} = 0.
    !> Occupancies are 0/1 per orbital (Pauli exclusion). The returned `factor`
    !> includes the parity sign resulting from reordering operators. If the
    !> excitation violates Pauli or leaves the active space, `iCNew` is set to 0.
    module subroutine ExciteConfiguration(this, iCNew, factor, creates, destroys, iC)
      class(T_ConfigList_E_AllActive_Fermionic), intent(in) :: this
      !> Output configuration index after excitation (0 if not possible).
      integer(I32), intent(out)        :: iCNew
      !> Output scalar factor including fermionic sign.
      real(R64), intent(out)           :: factor
      !> Orbital indices where particles are created (a^\dagger operators).
      integer(I32), intent(in), contiguous :: creates(:)
      !> Orbital indices where particles are annihilated (a operators).
      integer(I32), intent(in), contiguous :: destroys(:)
      !> Input configuration index to apply the excitation to.
      integer(I32), intent(in)         :: iC
    end subroutine

    !> Perform fermionic statistics-specific runtime setup. Typical tasks include
    !> allocating caches and finalizing mapping tables.
    module subroutine Setup(this)
      class(T_ConfigList_E_AllActive_Fermionic), intent(inout) :: this
    end subroutine

    !> Perform fermionic statistics-specific fabrication from input, preparing
    !> mapping structures between configuration indices and occupation encodings.
    module subroutine Fabricate(this)
      class(T_ConfigList_E_AllActive_Fermionic), intent(inout) :: this
    end subroutine
  end interface

end module
