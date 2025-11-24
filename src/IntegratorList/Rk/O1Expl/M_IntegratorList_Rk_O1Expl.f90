! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Runge–Kutta order 1 explicit (Forward Euler) backend.
!>
!> Concrete `T_IntegratorList_E_Rk` implementing the first-order
!> explicit Euler method. Useful as a simple baseline or when very small step
!> sizes are acceptable.
module M_IntegratorList_Rk_O1Expl
  use M_Utils_Types
  use M_IntegratorList, only: T_IntegratorList_E
  use M_IntegratorList_Rk, only: T_IntegratorList_E_Rk

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocate a first-order explicit RK integrator element.
    module subroutine IntegratorList_Rk_O1Expl_Allocate(e, path)
      class(T_IntegratorList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> First-order explicit RK (Euler) integrator element.
  type, extends(T_IntegratorList_E_Rk) :: T_IntegratorList_E_Rk_O1Expl
  contains
    !> Read variant-specific parameters from configuration.
    procedure :: FabricateLevel2 => Fabricate
    !> Allocate per-variant working memory.
    procedure :: SetupLevel2 => Setup
    !> Advance the state using the forward Euler step.
    procedure :: Integrate
  end type

  interface
    !> Fabricate the Euler integrator from configuration.
    module subroutine Fabricate(this)
      !> The Euler integrator instance to initialize
      class(T_IntegratorList_E_Rk_O1Expl), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Perform post-fabrication setup (allocate buffers, constants).
    module subroutine Setup(this)
      !> The Euler integrator instance to set up
      class(T_IntegratorList_E_Rk_O1Expl), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Integrate from `t0` to `t1` using the forward Euler method.
    !> Uses the `TimeDerivative` callback for RHS evaluations.
    module subroutine Integrate(this, state, t0, t1)
      !> The Euler integrator instance
      class(T_IntegratorList_E_Rk_O1Expl), intent(inout) :: this
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
