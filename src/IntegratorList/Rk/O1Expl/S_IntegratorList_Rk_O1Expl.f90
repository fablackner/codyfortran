! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_IntegratorList_Rk_O1Expl) S_IntegratorList_Rk_O1Expl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine IntegratorList_Rk_O1Expl_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_IntegratorList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_IntegratorList_E_Rk_O1Expl :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_IntegratorList_E_Rk_O1Expl), intent(inout) :: this

    call Say_Fabricate(this % path//".o1Expl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Rk_O1Expl), intent(inout) :: this

    call Say_Setup(this % path//".o1Expl")

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Integrate(this, state, t0, t1)
    class(T_IntegratorList_E_Rk_O1Expl), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0       ! Starting time
    real(R64), intent(in) :: t1       ! Target ending time

    complex(R64), allocatable :: dState(:)
    complex(R64) :: dt

    allocate (dState, mold=state)

    ! Calculate time step
    dt = cmplx(t1 - t0, 0.0_R64, R64)

    ! Calculate derivative at current time
    call this % TimeDerivative(dState, state, t0)

    ! Update state using forward Euler method
    state(:) = state(:) + dState(:) * dt

    deallocate (dState)

  end subroutine

end submodule
