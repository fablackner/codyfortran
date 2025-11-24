! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Runge–Kutta order 4 explicit (classical RK4) backend.
!>
!> Concrete `T_IntegratorList_E_Rk` implementing the classical
!> fourth-order explicit RK method, a widely used balance between accuracy and
!> cost for smooth problems.
module M_IntegratorList_Rk_O4Expl
  use M_Utils_Types
  use M_IntegratorList, only: T_IntegratorList_E
  use M_IntegratorList_Rk, only: T_IntegratorList_E_Rk

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocate a fourth-order explicit RK integrator element.
    module subroutine IntegratorList_Rk_O4Expl_Allocate(e, path)
      class(T_IntegratorList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> Fourth-order explicit RK (RK4) integrator element.
  type, extends(T_IntegratorList_E_Rk) :: T_IntegratorList_E_Rk_O4Expl
  contains
    !> Read variant-specific parameters from configuration.
    procedure :: FabricateLevel2 => Fabricate
    !> Allocate per-variant working memory.
    procedure :: SetupLevel2 => Setup
    !> Advance the state using the classical RK4 step.
    procedure :: Integrate
  end type

  interface
    !> Fabricate the RK4 integrator from configuration.
    module subroutine Fabricate(this)
      !> The RK4 integrator instance to initialize
      class(T_IntegratorList_E_Rk_O4Expl), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Perform post-fabrication setup (allocate buffers, constants).
    module subroutine Setup(this)
      !> The RK4 integrator instance to set up
      class(T_IntegratorList_E_Rk_O4Expl), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Integrate from `t0` to `t1` using the classical RK4 method.
    !> Uses the `TimeDerivative` callback for RHS evaluations.
    module subroutine Integrate(this, state, t0, t1)
      !> The RK4 integrator instance
      class(T_IntegratorList_E_Rk_O4Expl), intent(inout) :: this
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
