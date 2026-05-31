! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_IntegratorList_CrankNicolson.f90
!> @brief Crank–Nicolson integrator implementation.
!>
!> @details
!> Implements the second-order, A-stable Crank–Nicolson scheme:
!>
!>     ψ_{n+1} = ψ_n + (Δt/2) [f(t_n, ψ_n) + f(t_{n+1}, ψ_{n+1})]
!>
!> Because the scheme is implicit, the nonlinear equation is solved via
!> **fixed-point iteration** (Picard iteration). Starting from an explicit
!> Euler guess, we iterate until the residual drops below `tolerance` or
!> `maxIter` is reached.
!>
!> ## Properties
!>
!> - **Order**: 2
!> - **Stability**: A-stable (unconditionally stable for linear problems)
!> - **Symplecticity**: Norm-preserving for skew-Hermitian generators
!> - **Cost**: Requires iterative solve; 2 RHS evaluations per iteration
!>
!> ## Limitations
!>
!> The current implementation uses a simple fixed-point iteration which may
!> converge slowly for stiff or highly nonlinear problems. For such cases
!> consider the Expokit or GSL backends with adaptive stepping.
submodule(M_IntegratorList_Cn) S_IntegratorList_Cn

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine IntegratorList_Cn_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_IntegratorList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_IntegratorList_E_Cn :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_IntegratorList_E_Cn), intent(inout) :: this

    call Say_Fabricate(this % path//".cn")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Cn), intent(inout) :: this

    call Say_Setup(this % path//".cn")

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Advance the state using the Crank–Nicolson scheme.
  !>
  !> Solves the implicit equation via fixed-point iteration:
  !> 1. Initial guess: forward Euler (ψ_new = ψ + Δt·k1)
  !> 2. Iterate: ψ_new = ψ + (Δt/2)·(k1 + k2), where k2 = f(t1, ψ_new)
  !> 3. Stop when ||ψ_new - ψ_prev||_∞ < tolerance or maxIter reached.
  !>
  !> @warning Emits a console warning if iteration does not converge.
  module subroutine Integrate(this, state, t0, t1)

    class(T_IntegratorList_E_Cn), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1

    ! Work arrays for iterative solve
    complex(R64), allocatable :: k1(:), k2(:), stateNew(:), statePrev(:)
    real(R64) :: dt
    integer :: iter
    real(R64) :: error
    logical :: convergedQ

    ! Iteration parameters (could be made configurable via JSON)
    integer, parameter :: maxIter = 50
    real(R64), parameter :: tolerance = 1.0E-10_R64

    ! Allocate temporary arrays
    allocate (k1(size(state)))
    allocate (k2(size(state)))
    allocate (stateNew(size(state)))
    allocate (statePrev(size(state)))

    ! Calculate time step
    dt = t1 - t0

    ! Evaluate derivative at initial time: k1 = f(t0, state)
    call this % TimeDerivative(k1, state, t0)

    ! Initial guess using forward Euler
    stateNew = state + cmplx(dt, 0.0_R64, R64) * k1

    ! Fixed-point iteration for implicit part
    convergedQ = .false.

    do iter = 1, maxIter
      statePrev = stateNew

      ! Evaluate derivative at final time with current estimate: k2 = f(t1, stateNew)
      call this % TimeDerivative(k2, stateNew, t1)

      ! Crank–Nicolson update: ψ_{n+1} = ψ_n + (Δt/2)·(k1 + k2)
      stateNew = state + 0.5_R64 * cmplx(dt, 0.0_R64, R64) * (k1 + k2)

      ! Check convergence: ||stateNew - statePrev||_∞
      error = maxval(abs(stateNew - statePrev))

      if (error < tolerance) then
        convergedQ = .true.
        exit
      end if
    end do

    if (.not. convergedQ) then
      write (*, '(A,I0,A,ES10.3)') &
        "Warning: Crank-Nicolson iteration did not converge after ", maxIter, &
        " iterations. Residual = ", error
    end if

    ! Commit the solution
    state = stateNew

    deallocate (k1, k2, stateNew, statePrev)
  end subroutine

end submodule
