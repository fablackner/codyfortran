! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  module subroutine Integrate(this, state, t0, t1)
    class(T_IntegratorList_E_Rk_O4Expl), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0       ! Starting time
    real(R64), intent(in) :: t1       ! Target ending time

    complex(R64), allocatable :: stateTmp(:)
    complex(R64), allocatable :: increment(:)
    complex(R64), allocatable :: dState(:)
    real(R64) :: timePrime
    real(R64) :: dt, dt1, dt2, dt3, dt6

    allocate (dState, mold=state)
    allocate (increment, mold=state)
    allocate (stateTmp, mold=state)

    ! Calculate time step size
    dt = (t1 - t0)
    dt1 = dt
    dt2 = dt / 2.0_R64
    dt3 = dt / 3.0_R64
    dt6 = dt / 6.0_R64

    ! first step
    Call this % TimeDerivative(dState, state, t0)

    increment(:) = dState(:) * dt6
    stateTmp(:) = state(:) + dState(:) * dt2

    ! second step
    timePrime = t0 + dt * 0.5_R64
    Call this % TimeDerivative(dState, stateTmp, timePrime)

    increment(:) = increment(:) + dState(:) * dt3
    stateTmp(:) = state(:) + dState(:) * dt2

    ! third step
    Call this % TimeDerivative(dState, stateTmp, timePrime)

    increment(:) = increment(:) + dState(:) * dt3
    stateTmp(:) = state(:) + dState(:) * dt1

    ! fourth step
    timePrime = t0 + dt
    Call this % TimeDerivative(dState, stateTmp, timePrime)

    increment(:) = increment(:) + dState(:) * dt6
    state(:) = state(:) + increment(:)

    deallocate (dState)
    deallocate (increment)
    deallocate (stateTmp)

  end subroutine

end submodule
