! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Crank–Nicolson integrator backend.
!>
!> This module provides a concrete `T_IntegratorList_E` that advances the
!> state using the second-order, A-stable Crank–Nicolson scheme (the trapezoidal
!> rule / implicit midpoint for linear systems). It is time-reversible and
!> norm-preserving for skew-Hermitian generators.
module M_IntegratorList_Cn
  use M_Utils_Types
  use M_IntegratorList, only: T_IntegratorList_E

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocate a Crank–Nicolson integrator element.
    !>
    !> The resulting object will later be fabricated and set up using its
    !> type-bound procedures.
    module subroutine IntegratorList_Cn_Allocate(e, path)
      class(T_IntegratorList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> Crank–Nicolson integrator element.
  !>
  !> Implements `Fabricate`, `Setup`, and `Integrate` for the CN scheme.
  type, extends(T_IntegratorList_E) :: T_IntegratorList_E_Cn
  contains
    !> Read CN-specific parameters from configuration and capture dependencies.
    procedure :: Fabricate
    !> Allocate working memory and precompute constants.
    procedure :: Setup
    !> Advance the state from `t0` to `t1` using the Crank–Nicolson stepper.
    procedure :: Integrate
  end type

  interface
    !> Fabricate the CN integrator from configuration.
    module subroutine Fabricate(this)
      !> The CN integrator instance to initialize
      class(T_IntegratorList_E_Cn), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Perform post-fabrication setup (allocate buffers, constants).
    module subroutine Setup(this)
      !> The CN integrator instance to set up
      class(T_IntegratorList_E_Cn), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Integrate from `t0` to `t1` using the Crank–Nicolson scheme.
    !> Uses the `TimeDerivative` callback stored in the instance.
    module subroutine Integrate(this, state, t0, t1)
      !> The CN integrator instance
      class(T_IntegratorList_E_Cn), intent(inout) :: this
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
