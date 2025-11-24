! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Propagator_SplitStep_Order2) S_Propagator_SplitStep_Order2

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  subroutine Propagate(state, t0, t1)
    use M_IntegratorList
    use M_Propagator

    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in)    :: t0      ! Starting time
    real(R64), intent(in)    :: t1      ! Target ending time

    real(R64), parameter :: a1 = 0.5_R64
    real(R64), parameter :: b1 = 1.0_R64
    real(R64), parameter :: a2 = 0.5_R64

    ! Local variables
    real(R64)    :: dt, time1, time2

    ! Calculate time step (single step)
    dt = (t1 - t0)

    time1 = t0
    time2 = time1 + a1 * dt

    ! No need to pass TimeDerivative anymore as it was set in Setup
    call IntegratorList(1) % e % integrate(state, time1, time1 + a1 * dt)
    call IntegratorList(2) % e % integrate(state, time1, time1 + b1 * dt)
    call IntegratorList(1) % e % integrate(state, time2, time2 + a2 * dt)

  end subroutine

end submodule
