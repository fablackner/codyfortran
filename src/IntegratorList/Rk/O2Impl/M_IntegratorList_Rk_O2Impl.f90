! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Runge–Kutta order 2 implicit backend.
!>
!> Concrete `T_IntegratorList_E_Rk` implementing a second-order
!> implicit RK scheme (e.g., implicit midpoint or trapezoidal rule). Implicit
!> methods offer improved stability, especially for stiff problems.
module M_IntegratorList_Rk_O2Impl
  use M_Utils_Types
  use M_IntegratorList, only: T_IntegratorList_E
  use M_IntegratorList_Rk, only: T_IntegratorList_E_Rk

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocate a second-order implicit RK integrator element.
    module subroutine IntegratorList_Rk_O2Impl_Allocate(e, path)
      class(T_IntegratorList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> Second-order implicit RK integrator element.
  type, extends(T_IntegratorList_E_Rk) :: T_IntegratorList_E_Rk_O2Impl
  contains
    !> Read variant-specific parameters from configuration.
    procedure :: FabricateLevel2 => Fabricate
    !> Allocate per-variant working memory.
    procedure :: SetupLevel2 => Setup
    !> Advance the state using a second-order implicit RK step.
    procedure :: Integrate
  end type

  interface
    !> Fabricate the second-order implicit RK integrator from configuration.
    module subroutine Fabricate(this)
      !> The second-order implicit RK integrator instance to initialize
      class(T_IntegratorList_E_Rk_O2Impl), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Perform post-fabrication setup (allocate buffers, constants).
    module subroutine Setup(this)
      !> The second-order implicit RK integrator instance to set up
      class(T_IntegratorList_E_Rk_O2Impl), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Integrate from `t0` to `t1` using a second-order implicit RK method.
    !> Uses the `TimeDerivative` callback for RHS evaluations.
    module subroutine Integrate(this, state, t0, t1)
      !> The second-order implicit RK integrator instance
      class(T_IntegratorList_E_Rk_O2Impl), intent(inout) :: this
      !> Quantum state to propagate (modified in-place)
      complex(R64), intent(inout), contiguous :: state(:)
      !> Initial time
      real(R64), intent(in) :: t0
      !> Final time
      real(R64), intent(in) :: t1
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
