! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Sb) S_Method_Sb

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Sb_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Method

    call Say_Fabricate("method.sb")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Method_Setup => Setup
    Method_GetEnergy => GetEnergy
    Method_TimeDerivative => TimeDerivative

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_OrbsInit
    use M_Method

    call Say_Setup("method.sb")

    allocate (Method_state(1:Grid_nPoints))

    ! Initialize the single body orbital
    call OrbsInit_InitializeOrb(Method_state, 1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function GetEnergy(time) result(res)
    use M_Utils_Constants
    use M_Grid
    use M_Method

    real(R64)                :: res
    real(R64), intent(in)    :: time

    complex(R64), allocatable :: stateTmp(:)
    integer(I32) :: nG

    nG = Grid_nPoints

    ! Allocate temporary orbital for H|ψ>
    allocate (stateTmp(nG))

    call TimeDerivative(stateTmp, Method_state, time)
    stateTmp = IU * stateTmp  ! H|ψ> = i d|ψ>/dt

    ! Calculate energy as expectation value <ψ|H|ψ>
    res = real(Grid_InnerProduct(Method_state, stateTmp), kind=R64)

    deallocate (stateTmp)
  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivative(dState, state, time)
    use M_Utils_Constants
    use M_SysKinetic
    use M_SysPotential

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), allocatable :: externalPotential(:)
    complex(R64), allocatable :: stateTmp(:)

    dState(:) = 0.0_R64

    allocate (stateTmp, mold=state)

    ! For single body, only apply kinetic and potential operators
    ! No meanfield or exchange operators needed
    call SysKinetic_MultiplyWithKineticOp(dState, state, time)
    call SysPotential_FillExternalPotential(externalPotential, time)
    call SysPotential_MultiplyWithExternalPotential(stateTmp, externalPotential, state)

    dState(:) = dState(:) + stateTmp(:)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

end submodule
