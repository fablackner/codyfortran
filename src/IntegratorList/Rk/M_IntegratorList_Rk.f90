! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Runge–Kutta integrator family (base type).
!>
!> This module defines an abstract extension of `T_IntegratorList_E`
!> specialized for Runge–Kutta schemes. Concrete variants (order 1 explicit,
!> order 2 explicit/implicit, order 4 explicit, …) live in submodules and
!> implement the algorithm-specific details while reusing the shared lifecycle.
module M_IntegratorList_Rk
  use M_Utils_Types
  use M_IntegratorList, only: T_IntegratorList_E

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocate a Runge–Kutta integrator element (concrete subtype).
    module subroutine IntegratorList_Rk_Allocate(e, path)
      class(T_IntegratorList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> Abstract base for Runge–Kutta integrators.
  type, abstract, extends(T_IntegratorList_E) :: T_IntegratorList_E_Rk
    !> Order of the RK method (controls accuracy/cost). Common values:
    !> 1 (Euler), 2 (midpoint/Heun/implicit midpoint), 4 (classical RK4).
    integer(I32) :: order
    !> Whether the scheme is implicit (.true.) or explicit (.false.).
    logical :: implicitQ
  contains
    !> Read generic RK parameters from configuration and delegate specifics.
    procedure :: Setup
    !> Perform post-fabrication setup and delegate to `SetupLevel2`.
    procedure :: Fabricate
    !> Implementation-specific setup (dimensions, caches, etc.).
    procedure(I_SetupLevel2), deferred :: SetupLevel2
    !> Implementation-specific fabrication (read per-variant settings).
    procedure(I_FabricateLevel2), deferred :: FabricateLevel2
  end type

  interface
    !> Fabricate a generic RK integrator and delegate variant-specific logic.
    module subroutine Fabricate(this)
      !> The RK integrator instance to initialize
      class(T_IntegratorList_E_Rk), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Perform post-fabrication setup and delegate to `SetupLevel2`.
    module subroutine Setup(this)
      !> The RK integrator instance to set up
      class(T_IntegratorList_E_Rk), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    !> Variant-specific setup hook for RK integrators.
    subroutine I_SetupLevel2(this)
      import :: T_IntegratorList_E_Rk
      class(T_IntegratorList_E_Rk), intent(inout) :: this
    end subroutine

    !> Variant-specific fabrication hook for RK integrators.
    subroutine I_FabricateLevel2(this)
      import :: T_IntegratorList_E_Rk
      class(T_IntegratorList_E_Rk), intent(inout) :: this
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
