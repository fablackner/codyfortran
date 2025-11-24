! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_IntegratorList_Rk_O2Impl) S_IntegratorList_Rk_O2Impl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine IntegratorList_Rk_O2Impl_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_IntegratorList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_IntegratorList_E_Rk_O2Impl :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_IntegratorList_E_Rk_O2Impl), intent(inout) :: this

    call Say_Fabricate(this % path//".o2Impl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Rk_O2Impl), intent(inout) :: this

    call Say_Setup(this % path//".o2Impl")

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Integrate(this, state, t0, t1)
    class(T_IntegratorList_E_Rk_O2Impl), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0       ! Starting time
    real(R64), intent(in) :: t1       ! Target ending time

    complex(R64), allocatable :: stateNew(:)
    complex(R64), allocatable :: stateMid(:)
    complex(R64), allocatable :: dState(:)
    real(R64) :: midTime
    real(R64) :: dt
    integer :: j
    integer, parameter :: maxIter = 10
    real(R64), parameter :: tolerance = 1.0e-10_R64
    logical :: convergedQ

    allocate (stateNew, mold=state)
    allocate (stateMid, mold=state)
    allocate (dState, mold=state)

    ! Calculate time step
    dt = (t1 - t0)

    ! Initialize with explicit Euler as first guess
    call this % TimeDerivative(dState, state, t0)
    stateNew(:) = state(:) + dt * dState(:)

    ! Iterative solution for implicit midpoint method
    convergedQ = .false.
    midTime = t0 + 0.5_R64 * dt

    do j = 1, maxIter
      ! Calculate midpoint state
      stateMid(:) = 0.5_R64 * (state(:) + stateNew(:))

      ! Calculate derivative at midpoint
      call this % TimeDerivative(dState, stateMid, midTime)

      ! Update new state
      stateNew(:) = state(:) + dt * dState(:)

      ! Check for convergence
      if (maxval(abs(stateNew(:) - stateMid(:) * 2.0_R64 + state(:))) < tolerance) then
        convergedQ = .true.
        exit
      end if
    end do

    if (.not. convergedQ) then
      print *, "Warning: Implicit Runge-Kutta did not converge at t =", t0
    end if

    ! Update state
    state(:) = stateNew(:)

    deallocate (stateNew)
    deallocate (stateMid)
    deallocate (dState)

  end subroutine

end submodule
