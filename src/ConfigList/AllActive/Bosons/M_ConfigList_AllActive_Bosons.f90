! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> All-active bosonic configuration implementation.
!>
!> Concrete implementation of the all-active configuration element for bosons.
!> Provides the excitation application logic and statistics-specific setup and
!> fabrication steps. Bosonic operators obey commutation relations and do not
!> introduce fermionic sign changes; occupancy per orbital is unbounded within
!> the active space.
module M_ConfigList_AllActive_Bosonic
  use M_Utils_Types
  use M_ConfigList, only: T_ConfigList_E
  use M_ConfigList_AllActive, only: T_ConfigList_E_AllActive

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocate a bosonic all-active configuration element from input at `path`.
    !>
    !> Returns the concrete element as a polymorphic
    !> `class(T_ConfigList_E)` to be stored in the global list.
    module subroutine ConfigList_E_AllActive_Bosonic_Allocate(e, path)
      class(T_ConfigList_E), allocatable, intent(out) :: e
      !> Path into the input document for this element's configuration.
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> Concrete all-active bosonic configuration element.
  !>
  !> Implements bosonic excitation logic (commuting operators) and provides the
  !> statistics-specific level-2 setup/fabrication hooks used by the base class.
  type, extends(T_ConfigList_E_AllActive) :: T_ConfigList_E_AllActive_Bosonic

  contains

    !> Apply creation and annihilation operators using bosonic commutation rules.
    procedure :: ExciteConfiguration

    !> Bosonic level-2 setup invoked by the all-active base `Setup`.
    procedure :: SetupLevel2 => Setup

    !> Bosonic level-2 fabrication invoked by the all-active base `Fabricate`.
    procedure :: FabricateLevel2 => Fabricate

  end type

  interface
    !> Apply bosonic excitation operators to transform a configuration index.
    !>
    !> Bosonic operators obey the commutation relations
    !> [a_i, a^\dagger_j] = \delta_{ij}, and there is no fermionic sign. The
    !> returned `factor` captures any scalar prefactors implied by the ladder
    !> operators for the given occupation pattern. If the excitation would leave
    !> the active space or is invalid, `iCNew` is set to 0.
    module subroutine ExciteConfiguration(this, iCNew, factor, creates, destroys, iC)
      class(T_ConfigList_E_AllActive_Bosonic), intent(in) :: this
      !> Output configuration index after excitation (0 if not possible).
      integer(I32), intent(out)        :: iCNew
      !> Output scalar factor (no fermionic sign for bosons).
      real(R64), intent(out)           :: factor
      !> Orbital indices where particles are created (a^\dagger operators).
      integer(I32), intent(in), contiguous :: creates(:)
      !> Orbital indices where particles are annihilated (a operators).
      integer(I32), intent(in), contiguous :: destroys(:)
      !> Input configuration index to apply the excitation to.
      integer(I32), intent(in)         :: iC
    end subroutine

    !> Perform bosonic statistics-specific runtime setup. Typical tasks include
    !> allocating caches and finalizing mapping tables.
    module subroutine Setup(this)
      class(T_ConfigList_E_AllActive_Bosonic), intent(inout) :: this
    end subroutine

    !> Perform bosonic statistics-specific fabrication from input, preparing
    !> mapping structures between configuration indices and occupation encodings.
    module subroutine Fabricate(this)
      class(T_ConfigList_E_AllActive_Bosonic), intent(inout) :: this
    end subroutine
  end interface

end module
