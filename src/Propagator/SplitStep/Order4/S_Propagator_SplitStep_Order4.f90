! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Propagator_SplitStep_Order4) S_Propagator_SplitStep_Order4

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  subroutine Propagate(state, t0, t1)
    use M_IntegratorList
    use M_Propagator

    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in)    :: t0      ! Starting time
    real(R64), intent(in)    :: t1      ! Target ending time

    real(R64), parameter    :: o1 = 1.0_R64 / (2.0_R64 - 2.0_R64**(1.0_R64 / 3.0_R64))
    real(R64), parameter    :: o2 = 1.0_R64 - 2.0_R64 * o1
    real(R64), parameter    :: a1 = 0.5_R64 * o1
    real(R64), parameter    :: b1 = o1
    real(R64), parameter    :: a2 = 0.5_R64 * (o1 + o2)
    real(R64), parameter    :: b2 = o2
    real(R64), parameter    :: a3 = 0.5_R64 * (o1 + o2)
    real(R64), parameter    :: b3 = o1
    real(R64), parameter    :: a4 = 0.5_R64 * o1

    real(R64)    :: dt, time1, time2, time3, time4, time5, time6

    ! Calculate time step (single step)
    dt = (t1 - t0)

    time1 = t0
    time2 = time1 + a1 * dt
    time3 = time2 + a2 * dt
    time4 = time3 + b2 * dt
    time5 = time4 + a3 * dt
    time6 = time5 + b3 * dt

    ! No need to pass TimeDerivative anymore as it was set in Setup
    call IntegratorList(1) % e % integrate(state, time1, time1 + a1 * dt)
    call IntegratorList(2) % e % integrate(state, time1, time1 + b1 * dt)
    call IntegratorList(1) % e % integrate(state, time2, time2 + a2 * dt)
    call IntegratorList(2) % e % integrate(state, time3, time3 + b2 * dt)
    call IntegratorList(1) % e % integrate(state, time4, time4 + a3 * dt)
    call IntegratorList(2) % e % integrate(state, time5, time5 + b3 * dt)
    call IntegratorList(1) % e % integrate(state, time6, time6 + a4 * dt)

  end subroutine

end submodule
