! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Exponential integrators for u' = L u + N(t, u) without operator splitting.
!>
!> L is the part of the generator that the method propagates exactly
!> (`Method_PropagateExactPart`, e.g. -i h0 on Ylm grids with
!> method.mb.orbBased.exactFieldFreeQ) and whose phi functions it provides
!> (`Method_ApplyExactPartPhi`); N is the sum of the TimeDerivative callbacks
!> of IntegratorList(1) and IntegratorList(2) (their Integrate procedures are
!> not used). Schemes (`"propagator.exponential.<scheme>"`):
!> - `krogstad4`: exponential Runge–Kutta of order 4 (one step, 4 N per step)
!> - `adams`: exponential Adams predictor–corrector (k steps, 2 N per step)
!> - `adamsLawson`: Adams predictor–corrector in the interaction picture of L
!>   (Lawson multistep, 2 N per step)
!>
!> This module also provides the building blocks shared by the schemes.
module M_Propagator_Exponential
  use M_Utils_Types

  implicit none

  interface
    !> Reads `"propagator.exponential"`, checks that the method provides the
    !> exact part and its phi functions in real time, and dispatches to the
    !> selected scheme.
    module subroutine Propagator_Exponential_Fabricate
    end subroutine

    !> Krogstad step (J. Comput. Phys. 203, 72 (2005)), order 4, from t0 to t1
    !> in place; also the startup step of the multistep schemes.
    module subroutine Propagator_Exponential_KrogstadStep(state, t0, t1)
      !> State, propagated in place.
      complex(R64), intent(inout), contiguous :: state(:)
      !> Start time.
      real(R64), intent(in) :: t0
      !> End time.
      real(R64), intent(in) :: t1
    end subroutine

    !> dState = N(time, state): the sum of the TimeDerivative callbacks of
    !> IntegratorList(1) and IntegratorList(2).
    module subroutine Propagator_Exponential_NonExactDerivative(dState, state, time)
      !> Output: derivative of everything except the exactly propagated part.
      complex(R64), intent(out), contiguous :: dState(:)
      !> Input: state.
      complex(R64), intent(in), contiguous :: state(:)
      !> Time.
      real(R64), intent(in) :: time
    end subroutine

    !> Matrix C that maps the values of a polynomial at the nodes θ_j to its
    !> scaled Taylor coefficients at θ = 0: p(θ) = sum_m (C values)_m θ^m / m!.
    module function Propagator_Exponential_TaylorCoefficients(nodes) result(res)
      !> Distinct interpolation nodes (in units of the step).
      real(R64), intent(in) :: nodes(:)
      !> size(nodes) x size(nodes); row m + 1 gives the m-th coefficient.
      real(R64), allocatable :: res(:, :)
    end function

    !> Inserts g as the newest column of a step history (newest in column 1),
    !> dropping the oldest once all columns are in use.
    module subroutine Propagator_Exponential_PushHistory(history, nHistory, g)
      !> History (state size x number of steps).
      complex(R64), intent(inout), contiguous :: history(:, :)
      !> Number of valid columns, updated.
      integer(I32), intent(inout) :: nHistory
      !> New value.
      complex(R64), intent(in), contiguous :: g(:)
    end subroutine

    !> Whether the step t0 -> t1 continues the previous step tLast - hLast ->
    !> tLast with the same step size (otherwise a multistep history is stale).
    module function Propagator_Exponential_ContinuesQ(t0, t1, tLast, hLast) result(res)
      !> Start and end time of the new step.
      real(R64), intent(in) :: t0, t1
      !> End time and step size of the previous step.
      real(R64), intent(in) :: tLast, hLast
      logical :: res
    end function
  end interface

end module
