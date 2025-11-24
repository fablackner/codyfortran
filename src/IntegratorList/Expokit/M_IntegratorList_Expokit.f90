! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Expokit-based Krylov exponential integrator backend.
!>
!> This module provides a concrete `T_IntegratorList_E` that advances the
!> state via actions of the matrix exponential using a Krylov (Arnoldi/Lanczos)
!> approximation, leveraging ideas from Expokit. It is well-suited for large,
!> sparse generators where forming the full exponential is infeasible.
module M_IntegratorList_Expokit
  use M_Utils_Types
  use M_IntegratorList, only: T_IntegratorList_E

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocate an Expokit-style Krylov integrator element.
    module subroutine IntegratorList_Expokit_Allocate(e, path)
      class(T_IntegratorList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> Expokit-style Krylov integrator element.
  type, extends(T_IntegratorList_E) :: T_IntegratorList_E_Expokit
    !> Krylov subspace dimension for the Arnoldi/Lanczos process.
    integer(I32) :: krylov_dim = 30
    !> Error tolerance for the matrix exponential approximation.
    real(R64) :: tolerance = 1.0e-7_R64
    !> Maximum number of internal steps in one `Integrate` call.
    integer(I32) :: max_steps = 1000
  contains
    !> Read parameters (Krylov dimension, tolerances, limits) from config.
    procedure :: Fabricate
    !> Allocate working memory and precompute workspace sizes.
    procedure :: Setup
    !> Advance the state using a Krylov approximation to exp(Δt·A)·v.
    procedure :: Integrate
  end type

  interface
    !> Fabricate the Krylov/Expokit integrator from configuration.
    module subroutine Fabricate(this)
      !> The Expokit-style integrator instance to initialize
      class(T_IntegratorList_E_Expokit), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Perform post-fabrication setup (allocate buffers, constants).
    module subroutine Setup(this)
      !> The Expokit-style integrator instance to set up
      class(T_IntegratorList_E_Expokit), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Integrate from `t0` to `t1` using Krylov-based exponential actions.
    !> Uses the `TimeDerivative` callback for matrix-vector products.
    module subroutine Integrate(this, state, t0, t1)
      !> The Expokit-style integrator instance
      class(T_IntegratorList_E_Expokit), intent(inout) :: this
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
