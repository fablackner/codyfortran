! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Implementation of the fourth-order split-operator propagator.
!>
!> Applies the Yoshida/Forest–Ruth symmetric composition:
!>   S₄(Δt) = S₂(γ₁Δt) · S₂(γ₂Δt) · S₂(γ₁Δt)
!>
!> where S₂ is the 2nd-order Strang splitting and the coefficients are:
!>   γ₁ = 1/(2 − 2^(1/3)) ≈ 1.3512
!>   γ₂ = 1 − 2γ₁ ≈ −1.7024  (negative "backward" step)
!>
!> Properties:
!> - Time-reversible (symmetric)
!> - Local error: O(Δt⁵)
!> - Global error: O(Δt⁴)
!> - 7 operator applications per step (vs. 3 for order-2)
!> - Requires 2 integrators in IntegratorList: (1) for Â, (2) for B̂
!>
!> Reference: Yoshida, Phys. Lett. A 150, 262 (1990)
submodule(M_Propagator_SplitStep_Order4) S_Propagator_SplitStep_Order4

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Binds the Propagator facade to the order-4 split-step implementation.
  module subroutine Propagator_SplitStep_Order4_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Propagator

    call Say_Fabricate("propagator.splitStep.order4")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Propagator_Propagate => Propagate

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Applies the 7-stage Yoshida symmetric composition.
!>
!> Sequence: A(a₁) B(b₁) A(a₂) B(b₂) A(a₃) B(b₃) A(a₄)
!>
!> The coefficients satisfy:
!>   Σaᵢ = Σbᵢ = 1  (consistency)
!>   Symmetry: a₁=a₄, a₂=a₃, b₁=b₃
!>   Fourth-order conditions via Yoshida's construction
!>
!> Arguments:
!>   state  - Complex state vector, modified in-place
!>   t0     - Initial time
!>   t1     - Final time
  subroutine Propagate(state, t0, t1)
    use M_IntegratorList
    use M_Propagator

    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in)    :: t0
    real(R64), intent(in)    :: t1

    ! Yoshida 4th-order coefficients: γ₁ = 1/(2 − 2^(1/3)), γ₂ = 1 − 2γ₁
    real(R64), parameter    :: o1 = 1.0_R64 / (2.0_R64 - 2.0_R64**(1.0_R64 / 3.0_R64))
    real(R64), parameter    :: o2 = 1.0_R64 - 2.0_R64 * o1

    ! Derived splitting coefficients: ABAB...A pattern (7 stages)
    real(R64), parameter    :: a1 = 0.5_R64 * o1           ! ≈ 0.6756
    real(R64), parameter    :: b1 = o1                     ! ≈ 1.3512
    real(R64), parameter    :: a2 = 0.5_R64 * (o1 + o2)    ! ≈ −0.1756
    real(R64), parameter    :: b2 = o2                     ! ≈ −1.7024
    real(R64), parameter    :: a3 = 0.5_R64 * (o1 + o2)    ! ≈ −0.1756 (symmetric)
    real(R64), parameter    :: b3 = o1                     ! ≈ 1.3512 (symmetric)
    real(R64), parameter    :: a4 = 0.5_R64 * o1           ! ≈ 0.6756 (symmetric)

    real(R64)    :: dt, time1, time2, time3, time4, time5, time6

    dt = (t1 - t0)

    ! Cumulative time markers for time-dependent operators
    time1 = t0
    time2 = time1 + a1 * dt
    time3 = time2 + a2 * dt
    time4 = time3 + b2 * dt
    time5 = time4 + a3 * dt
    time6 = time5 + b3 * dt

    ! 7-stage symmetric sequence: A B A B A B A
    call IntegratorList(1) % e % Integrate(state, time1, time1 + a1 * dt)
    call IntegratorList(2) % e % Integrate(state, time1, time1 + b1 * dt)
    call IntegratorList(1) % e % Integrate(state, time2, time2 + a2 * dt)
    call IntegratorList(2) % e % Integrate(state, time3, time3 + b2 * dt)
    call IntegratorList(1) % e % Integrate(state, time4, time4 + a3 * dt)
    call IntegratorList(2) % e % Integrate(state, time5, time5 + b3 * dt)
    call IntegratorList(1) % e % Integrate(state, time6, time6 + a4 * dt)

  end subroutine

end submodule
