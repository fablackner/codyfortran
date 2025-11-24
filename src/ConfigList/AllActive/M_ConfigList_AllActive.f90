! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> All-active configuration base and factory.
!>
!> This module defines the abstract extension of the core configuration element
!> tailored for the "all-active" truncation scheme (all orbitals in the active
!> space can participate up to a fixed excitation rank). It also declares the
!> allocation factory used by higher layers to create concrete bosonic or
!> fermionic implementations and return them via the common base type.
module M_ConfigList_AllActive
  use M_Utils_Types
  use M_ConfigList, only: T_ConfigList_E

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocate a concrete all-active configuration element.
    !>
    !> The implementation inspects the input located at `path` and allocates a
    !> concrete `T_ConfigList_E_AllActive_*` (bosonic/fermionic),
    !> returning it as a polymorphic `T_ConfigList_E`.
    module subroutine ConfigList_AllActive_Allocate(e, path)
      class(T_ConfigList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> Abstract base class for all-active configuration elements.
  !>
  !> Adds the maximum excitation rank and provides concrete `Setup`/`Fabricate`
  !> bindings that orchestrate common all-active logic and then delegate to the
  !> statistics-specific level-2 procedures.
  type, abstract, extends(T_ConfigList_E) :: T_ConfigList_E_AllActive

    !> Maximum excitation rank supported (e.g., 1 for singles, 2 for doubles).
    integer(I32)     :: nExcitations

  contains
    !> Orchestrate runtime setup for all-active elements and call `SetupLevel2`.
    procedure :: Setup
    !> Orchestrate fabrication from input for all-active elements and call `FabricateLevel2`.
    procedure :: Fabricate
    !> Statistics-specific second-stage setup performed by concrete subclasses.
    procedure(I_SetupLevel2), deferred :: SetupLevel2
    !> Statistics-specific second-stage fabrication performed by concrete subclasses.
    procedure(I_FabricateLevel2), deferred :: FabricateLevel2
  end type

  interface
    !> Entry point implemented in this module that performs all-active common
    !> setup and then calls the subclass' `SetupLevel2`.
    module subroutine Setup(this)
      class(T_ConfigList_E_AllActive), intent(inout) :: this
    end subroutine

    !> Entry point implemented in this module that performs all-active common
    !> fabrication and then calls the subclass' `FabricateLevel2`.
    module subroutine Fabricate(this)
      class(T_ConfigList_E_AllActive), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    !> Statistics-specific second-stage setup. Called by `Setup` after the
    !> all-active common work has been performed.
    subroutine I_SetupLevel2(this)
      import :: T_ConfigList_E_AllActive
      class(T_ConfigList_E_AllActive), intent(inout) :: this
    end subroutine

    !> Statistics-specific second-stage fabrication. Called by `Fabricate` after
    !> the all-active common work has been performed.
    subroutine I_FabricateLevel2(this)
      import :: T_ConfigList_E_AllActive
      class(T_ConfigList_E_AllActive), intent(inout) :: this
    end subroutine
  end interface

end module
