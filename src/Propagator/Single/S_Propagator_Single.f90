! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Implementation of the single-RHS propagator backend.
!>
!> Provides monolithic time integration by delegating directly to the first
!> integrator in `IntegratorList`. No operator splitting is applied.
!>
!> Requirements:
!> - IntegratorList must be fabricated and set up before propagation.
!> - IntegratorList(1) must have its TimeDerivative callback configured.
submodule(M_Propagator_Single) S_Propagator_Single

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Binds the Propagator facade pointers to this backend.
!>
!> After calling this routine:
!> - `Propagator_Propagate` → Propagate (below)
!> - `Propagator_Setup` remains at its default (no-op or prior binding)
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
!> Advances the quantum state from t0 to t1 via IntegratorList(1).
!>
!> The integrator (RK4, SIL, Expokit, etc.) handles step-size control and
!> calls Method_TimeDerivative internally to evaluate dΨ/dt = −iĤΨ.
!>
!> Arguments:
!>   state  - Complex state vector, modified in-place
!>   t0     - Initial time
!>   t1     - Final time (may be < t0 for backward propagation)
  subroutine Propagate(state, t0, t1)
    use M_IntegratorList
    use M_Propagator

    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in)    :: t0
    real(R64), intent(in)    :: t1

    call IntegratorList(1) % e % Integrate(state, t0, t1)

  end subroutine

end submodule
