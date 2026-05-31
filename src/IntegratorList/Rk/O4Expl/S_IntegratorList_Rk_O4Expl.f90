! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_IntegratorList_Rk_O4Expl.f90
!> @brief Classical fourth-order Runge–Kutta (RK4) implementation.
!>
!> @details
!> Implements the classical 4th-order explicit RK method with Butcher tableau:
!>
!> @verbatim
!>     0   |
!>    1/2  | 1/2
!>    1/2  |  0   1/2
!>     1   |  0    0    1
!>    -----+--------------------
!>         | 1/6  1/3  1/3  1/6
!> @endverbatim
!>
!> The update formula is:
!>
!>     k1 = f(t_n, y_n)
!>     k2 = f(t_n + Δt/2, y_n + Δt/2 · k1)
!>     k3 = f(t_n + Δt/2, y_n + Δt/2 · k2)
!>     k4 = f(t_n + Δt, y_n + Δt · k3)
!>     y_{n+1} = y_n + Δt/6 · (k1 + 2k2 + 2k3 + k4)
!>
!> ## Properties
!>
!> - **Order**: 4 (local truncation error O(Δt⁵))
!> - **Stability**: Conditionally stable; CFL-like constraint for oscillatory
!>   problems.
!> - **Cost**: 4 RHS evaluations per step.
!>
!> This is the workhorse explicit integrator for smooth, non-stiff problems.
submodule(M_IntegratorList_Rk_O4Expl) S_IntegratorList_Rk_O4Expl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine IntegratorList_Rk_O4Expl_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_IntegratorList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_IntegratorList_E_Rk_O4Expl :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_IntegratorList_E_Rk_O4Expl), intent(inout) :: this

    call Say_Fabricate(this % path//".o4Expl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Rk_O4Expl), intent(inout) :: this

    call Say_Setup(this % path//".o4Expl")

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Advance the state using the classical RK4 method.
  !>
  !> Performs four RHS evaluations and combines them with weights (1/6, 1/3,
  !> 1/3, 1/6) to achieve 4th-order accuracy.
  !>
  !> @param[inout] this   The RK4 integrator instance.
  !> @param[inout] state  Complex state vector (modified in-place).
  !> @param[in]    t0     Starting time.
  !> @param[in]    t1     Target ending time.
  module subroutine Integrate(this, state, t0, t1)
    class(T_IntegratorList_E_Rk_O4Expl), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1

    complex(R64), allocatable :: stateTmp(:)  ! Intermediate state
    complex(R64), allocatable :: increment(:) ! Accumulated weighted sum
    complex(R64), allocatable :: dState(:)    ! Current stage derivative
    real(R64) :: timePrime
    real(R64) :: dt, dt2, dt3, dt6

    allocate (dState, mold=state)
    allocate (increment, mold=state)
    allocate (stateTmp, mold=state)

    dt = t1 - t0
    dt2 = dt / 2.0_R64
    dt3 = dt / 3.0_R64
    dt6 = dt / 6.0_R64

    ! Stage 1: k1 = f(t0, state)
    call this % TimeDerivative(dState, state, t0)

    increment(:) = dState(:) * dt6            ! weight 1/6
    stateTmp(:) = state(:) + dState(:) * dt2 ! half-step for stage 2

    ! Stage 2: k2 = f(t0 + dt/2, state + dt/2 * k1)
    timePrime = t0 + dt2
    call this % TimeDerivative(dState, stateTmp, timePrime)

    increment(:) = increment(:) + dState(:) * dt3  ! weight 1/3
    stateTmp(:) = state(:) + dState(:) * dt2      ! half-step for stage 3

    ! Stage 3: k3 = f(t0 + dt/2, state + dt/2 * k2)
    call this % TimeDerivative(dState, stateTmp, timePrime)

    increment(:) = increment(:) + dState(:) * dt3  ! weight 1/3
    stateTmp(:) = state(:) + dState(:) * dt       ! full step for stage 4

    ! Stage 4: k4 = f(t0 + dt, state + dt * k3)
    timePrime = t0 + dt
    call this % TimeDerivative(dState, stateTmp, timePrime)

    increment(:) = increment(:) + dState(:) * dt6  ! weight 1/6

    ! Final update: y_{n+1} = y_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    state(:) = state(:) + increment(:)

    deallocate (dState, increment, stateTmp)

  end subroutine

end submodule
