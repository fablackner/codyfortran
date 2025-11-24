! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_IntegratorList_Rk_O2Expl) S_IntegratorList_Rk_O2Expl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine IntegratorList_Rk_O2Expl_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_IntegratorList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_IntegratorList_E_Rk_O2Expl :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_IntegratorList_E_Rk_O2Expl), intent(inout) :: this

    call Say_Fabricate(this % path//".o2Expl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Rk_O2Expl), intent(inout) :: this

    call Say_Setup(this % path//".o2Expl")

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Integrate(this, state, t0, t1)
    class(T_IntegratorList_E_Rk_O2Expl), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0       ! Starting time
    real(R64), intent(in) :: t1       ! Target ending time

    complex(R64), allocatable :: k1(:), k2(:), stateTmp(:)
    real(R64) :: dt

    allocate (k1, mold=state)
    allocate (k2, mold=state)
    allocate (stateTmp, mold=state)

    dt = (t1 - t0)

    ! First stage
    call this % TimeDerivative(k1, state, t0)

    ! Calculate intermediate state
    stateTmp(:) = state(:) + 0.5_R64 * dt * k1(:)

    ! Second stage
    call this % TimeDerivative(k2, stateTmp, t0 + 0.5_R64 * dt)

    ! Final update
    state(:) = state(:) + dt * k2(:)

    deallocate (k1)
    deallocate (k2)
    deallocate (stateTmp)

  end subroutine

end submodule
