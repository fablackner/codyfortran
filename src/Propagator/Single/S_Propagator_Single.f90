! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Propagator_Single) S_Propagator_Single

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Propagator_Single_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Propagator

    call Say_Fabricate("propagator.single")

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

    ! No need to pass TimeDerivative anymore as it was set in Setup
    call IntegratorList(1) % e % Integrate(state, t0, t1)

  end subroutine

end submodule
