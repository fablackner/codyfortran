! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_IntegratorList_Expokit) S_IntegratorList_Expokit

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine IntegratorList_Expokit_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_IntegratorList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_IntegratorList_E_Expokit :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_IntegratorList_E_Expokit), intent(inout) :: this

    call Say_Fabricate(this % path//".expokit")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    this % krylov_dim = Json_Get("krylov_dim", 30, path_=this % path)
    this % tolerance = Json_Get("tolerance", 1.0e-7_R64, path_=this % path)
    this % max_steps = Json_Get("max_steps", 1000, path_=this % path)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Expokit), intent(inout) :: this

    call Say_Setup(this % path//".expokit")

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Integrate(this, state, t0, t1)
    use M_Utils_ExpokitLib, only: ExpokitLib_IntegrateSym
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Expokit), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0       ! Starting time
    real(R64), intent(in) :: t1       ! Target ending time

    ! Local variables
    integer :: iflag, nSteps
    real(R64) :: t_step
    character(len=100) :: error_msg

    ! Calculate time step
    t_step = t1 - t0

    ! Set default number of steps for Expokit
    nSteps = this % max_steps

    ! Call the ExpokitLib wrapper for Hermitian matrices
    call ExpokitLib_IntegrateSym(state, t_step, ApplyMatrix, &
                                 this % krylov_dim, this % tolerance, nSteps, iflag)

    ! Check for errors
    if (iflag .ne. 0) then
      write (error_msg, '(A,I0)') "Expokit zhexpv failed with status: ", iflag
      error stop error_msg
    end if

  contains
    ! Matrix-vector product callback required by ExpokitLib
    ! Converts from the TimeDerivative interface (-i*H*x) to H*x for Expokit
    subroutine ApplyMatrix(dState, inState)
      complex(R64), intent(out), contiguous, target :: dState(:)
      complex(R64), intent(in), contiguous, target  :: inState(:)

      ! The TimeDerivative gives us -i*H*x, we need H*x for Expokit
      ! So we call TimeDerivative and then multiply by i
      call this % TimeDerivative(dState, inState, t0)
      dState = (0.0_R64, 1.0_R64) * dState  ! Multiply by i
    end subroutine
  end subroutine

end submodule
