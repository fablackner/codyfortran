! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Implementation of the second-order split-operator propagator.
!>
!> Applies the symmetric Strang splitting (a.k.a. leapfrog):
!>   exp(−iĤΔt) ≈ exp(−iÂ·Δt/2) · exp(−iB̂·Δt) · exp(−iÂ·Δt/2)
!>
!> Properties:
!> - Time-reversible (symmetric)
!> - Local error: O(Δt³)
!> - Global error: O(Δt²)
!> - Requires 2 integrators in IntegratorList: (1) for Â, (2) for B̂
!>
!> Typical use: Â = kinetic (diagonal in momentum), B̂ = potential (diagonal in position)
submodule(M_Propagator_SplitStep_Order2) S_Propagator_SplitStep_Order2

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Binds the Propagator facade to the order-2 split-step implementation.
  module subroutine Propagator_SplitStep_Order2_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Propagator

    call Say_Fabricate("propagator.splitStep.order2")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Propagator_Propagate => Propagate

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Applies the symmetric ABA splitting sequence.
!>
!> Sequence: A(½Δt) → B(Δt) → A(½Δt)
!>
!> Arguments:
!>   state  - Complex state vector, modified in-place
!>   t0     - Initial time
!>   t1     - Final time
!>
!> The time arguments passed to each integrator allow time-dependent operators.
!> IntegratorList(1) handles Â; IntegratorList(2) handles B̂.
  subroutine Propagate(state, t0, t1)
    use M_IntegratorList
    use M_Propagator

    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in)    :: t0
    real(R64), intent(in)    :: t1

    ! Strang splitting coefficients: A(½) B(1) A(½)
    real(R64), parameter :: a1 = 0.5_R64   ! First half-step of Â
    real(R64), parameter :: b1 = 1.0_R64   ! Full step of B̂
    real(R64), parameter :: a2 = 0.5_R64   ! Second half-step of Â

    real(R64)    :: dt, time1, time2

    dt = (t1 - t0)

    time1 = t0
    time2 = time1 + a1 * dt

    ! ABA sequence: apply Â for half-step, B̂ for full step, Â for half-step
    call IntegratorList(1) % e % Integrate(state, time1, time1 + a1 * dt)
    call IntegratorList(2) % e % Integrate(state, time1, time1 + b1 * dt)
    call IntegratorList(1) % e % Integrate(state, time2, time2 + a2 * dt)

  end subroutine

end submodule
