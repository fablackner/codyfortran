! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  module subroutine Integrate(this, state, t0, t1)

    class(T_IntegratorList_E_Cn), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0       ! Starting time
    real(R64), intent(in) :: t1       ! Target ending time

    complex(R64), allocatable :: k1(:), k2(:), stateNew(:), statePrev(:)
    real(R64) :: dt
    integer :: iter, maxIter
    real(R64) :: tolerance, error
    logical :: convergedQ

    ! Set parameters for iterative solver
    maxIter = 50
    tolerance = 1.0E-10_R64

    ! Allocate temporary arrays
    allocate (k1(size(state)))
    allocate (k2(size(state)))
    allocate (stateNew(size(state)))
    allocate (statePrev(size(state)))

    ! Calculate time step
    dt = t1 - t0

    ! Calculate first derivative at current time
    call this % TimeDerivative(k1, state, t0)

    ! Initial guess using forward Euler
    stateNew = state + cmplx(dt, 0.0_R64, R64) * k1

    ! Iterative solution for implicit part
    convergedQ = .false.

    do iter = 1, maxIter
      statePrev = stateNew

      ! Calculate second derivative at next time with current estimate
      call this % TimeDerivative(k2, stateNew, t1)

      ! Crank-Nicolson formula: y_{n+1} = y_n + (dt/2)*(k1 + k2)
      stateNew = state + 0.5_R64 * cmplx(dt, 0.0_R64, R64) * (k1 + k2)

      ! Check convergence
      error = maxval(abs(stateNew - statePrev))

      if (error < tolerance) then
        convergedQ = .true.
        exit
      end if
    end do

    if (.not. convergedQ) then
      write (*, *) "Warning: Crank-Nicolson iteration did not converge. Error =", error
    end if

    ! Update state with final solution
    state = stateNew

    ! Clean up
    deallocate (k1, k2, stateNew, statePrev)
  end subroutine

end submodule
