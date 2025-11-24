! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Short Iterative Lanczos (SIL) propagator backend.
!>
!> This module provides a concrete `T_IntegratorList_E` that advances the
!> state using a short Krylov subspace built via the Lanczos process. It
!> approximates the action of the exponential with low memory cost and is
!> effective for Hermitian/skew-Hermitian generators.
module M_IntegratorList_Sil
  use M_Utils_Types
  use M_IntegratorList, only: T_IntegratorList_E

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocate a Short Iterative Lanczos integrator element.
    module subroutine IntegratorList_Sil_Allocate(e, path)
      class(T_IntegratorList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> Short Iterative Lanczos integrator element.
  type, extends(T_IntegratorList_E) :: T_IntegratorList_E_Sil
  contains
    !> Read SIL-specific parameters from configuration and capture dependencies.
    procedure :: Fabricate
    !> Allocate working memory and precompute constants.
    procedure :: Setup
    !> Advance the state using a short Krylov/Lanczos exponential action.
    procedure :: Integrate
  end type

  interface
    !> Fabricate the SIL integrator from configuration.
    module subroutine Fabricate(this)
      !> The SIL integrator instance to initialize
      class(T_IntegratorList_E_Sil), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Perform post-fabrication setup (allocate buffers, constants).
    module subroutine Setup(this)
      !> The SIL integrator instance to set up
      class(T_IntegratorList_E_Sil), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Integrate from `t0` to `t1` using a short iterative Lanczos scheme.
    !> Uses the `TimeDerivative` callback for operator-vector actions.
    module subroutine Integrate(this, state, t0, t1)
      !> The SIL integrator instance
      class(T_IntegratorList_E_Sil), intent(inout) :: this
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
