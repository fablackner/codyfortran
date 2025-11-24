! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_IntegratorList_GslOdeiv2) S_IntegratorList_GslOdeiv2

  implicit none

  ! Parameters for GSL integration
  real(R64), parameter :: DEFAULT_ABS_ERROR = 1.0e-6_R64
  real(R64), parameter :: DEFAULT_REL_ERROR = 1.0e-10_R64

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine IntegratorList_GslOdeiv2_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_IntegratorList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_IntegratorList_E_GslOdeiv2 :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_IntegratorList_E_GslOdeiv2), intent(inout) :: this

    call Say_Fabricate(this % path//".gslOdeiv2")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    this % stepperType = Json_Get("stepperType", "rk4", path_=this % path)
    this % initializedQ = .false.  ! Will be initialized on first call to Integrate

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_GslOdeiv2), intent(inout) :: this

    call Say_Setup(this % path//".gslOdeiv2")

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Integrate(this, state, t0, t1)
    use M_Utils_Odeiv2GslLib

    class(T_IntegratorList_E_GslOdeiv2), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0       ! Starting time
    real(R64), intent(in) :: t1       ! Target ending time

    integer(I32) :: status
    integer :: nSteps

    ! Calculate dimension and step_size on first integration
    if (.not. this % initializedQ) then
      this % dimension = size(state)
      this % step_size = (t1 - t0)

      ! Initialize GSL with the stored TimeDerivative procedure
      call Odeiv2GslLib_CreateCtx(this % odeiv2Ctx, this % dimension, this % TimeDerivative, &
                                  this % step_size, DEFAULT_ABS_ERROR, DEFAULT_REL_ERROR, &
                                  trim(this % stepperType))

      this % initializedQ = .true.
    end if

    ! Set default number of steps
    nSteps = 1

    call Odeiv2GslLib_Integrate(status, state, t0, t1, nSteps, this % odeiv2Ctx)

  end subroutine

end submodule
