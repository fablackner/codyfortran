! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Adams–Lawson PECE step in the interaction picture v = E(tRef - t) u:
!>   P: vP      = v_n + h sum_j ab_j G_{n+1-j}                (k-step Adams–Bashforth)
!>   E: GP      = G(t_{n+1}, vP)
!>   C: v_{n+1} = v_n + h (am_0 GP + sum_j am_j G_{n+1-j})    (k-step Adams–Moulton)
!>   E: G_{n+1} = G(t_{n+1}, v_{n+1})
!> Every evaluation of G costs one N and two exact propagations; the last one
!> also returns u_{n+1} = E(t_{n+1} - tRef) v_{n+1}. The Adams weights are
!> the integrals over [0, 1] of the Lagrange polynomials of the nodes
!> (Propagator_Exponential_TaylorCoefficients). The history is built with
!> Krogstad steps on the first call and whenever a step does not continue the
!> previous one. If the caller changed the state between steps (absorber),
!> v_n is recomputed from it (one more exact propagation); the history is kept.
submodule(M_Propagator_Exponential_AdamsLawson) S_Propagator_Exponential_AdamsLawson

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Number of history values k (propagator.exponential.adamsLawson.nSteps).
  integer(I32) :: nSteps

  !> Adams–Bashforth weights (k) and Adams–Moulton weights (k + 1, new point first).
  real(R64), allocatable :: ab(:)
  real(R64), allocatable :: am(:)

  !> History G_n, G_{n-1}, ... (newest in column 1) and its number of valid columns.
  complex(R64), allocatable :: history(:, :)
  integer(I32) :: nHistory = 0

  !> Reference time of the interaction picture, end time and step size of the last step.
  real(R64) :: tRef, tLast, hLast

  !> State and interaction-picture state at tLast.
  complex(R64), allocatable :: uLast(:)
  complex(R64), allocatable :: vLast(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Reads `"propagator.exponential.adamsLawson.nSteps"` (1..6, default 4: order 5).
  module subroutine Propagator_Exponential_AdamsLawson_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Propagator
    use M_Propagator_Exponential

    integer(I32) :: j

    call Say_Fabricate("propagator.exponential.adamsLawson")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    nSteps = Json_Get("propagator.exponential.adamsLawson.nSteps", 4)
    if (nSteps < 1 .or. nSteps > 6) error stop "propagator.exponential.adamsLawson.nSteps must be in 1..6"

    ab = AdamsWeights([(real(1 - j, kind=R64), j=1, nSteps)])
    am = AdamsWeights([(real(2 - j, kind=R64), j=1, nSteps + 1)])
    nHistory = 0

    Propagator_Propagate => Propagate

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Weights of the quadrature ∫_0^1 p(θ) dθ for the interpolation polynomial
  !> p through the nodes: sum_m C_{mj} / (m + 1)!.
  function AdamsWeights(nodes) result(res)
    use M_Propagator_Exponential

    real(R64), intent(in) :: nodes(:)
    real(R64), allocatable :: res(:)

    real(R64), allocatable :: c(:, :)
    integer(I32) :: m

    c = Propagator_Exponential_TaylorCoefficients(nodes)
    allocate (res(size(nodes)))
    res(:) = 0.0_R64
    do m = 0, size(nodes) - 1
      res(:) = res(:) + c(m + 1, :) / gamma(real(m + 2, kind=R64))
    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Propagate(state, t0, t1)
    use M_Method
    use M_Propagator_Exponential

    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1

    complex(R64), allocatable :: v(:), vNew(:), g(:)
    real(R64) :: h
    integer(I32) :: j
    logical :: restartQ

    h = t1 - t0
    ! (re)allocate for the first call, a new nSteps or a new state size
    if (allocated(history)) then
      if (size(history, 1) .ne. size(state) .or. size(history, 2) .ne. nSteps) deallocate (history, uLast, vLast)
    end if
    if (.not. allocated(history)) then
      allocate (history(size(state), nSteps))
      allocate (uLast, vLast, mold=state)
      nHistory = 0
    end if
    allocate (vNew, g, mold=state)

    if (nHistory > 0) then
      if (.not. Propagator_Exponential_ContinuesQ(t0, t1, tLast, hLast)) nHistory = 0
    end if

    restartQ = nHistory .eq. 0
    if (restartQ) then
      ! G_n = N_n at the new reference time tRef = t0
      tRef = t0
      call Propagator_Exponential_NonExactDerivative(g, state, t0)
      call Propagator_Exponential_PushHistory(history, nHistory, g)
    end if

    !------------------------------------
    ! step: Krogstad until the history is complete, then PECE
    !------------------------------------

    if (nHistory < nSteps) then
      call Propagator_Exponential_KrogstadStep(state, t0, t1)
      vNew(:) = state(:)
      call Method_PropagateExactPart(vNew, t1, tRef)
    else
      ! interaction-picture state at t0 (vLast belongs to an older tRef after a restart)
      if (restartQ) then
        v = state
      else if (all(state .eq. uLast)) then
        v = vLast
      else
        v = state
        call Method_PropagateExactPart(v, t0, tRef)
      end if

      ! predictor
      vNew(:) = v(:)
      do j = 1, nSteps
        vNew(:) = vNew(:) + (h * ab(j)) * history(:, j)
      end do
      call EvaluateG(g, state, vNew, t1)

      ! corrector
      vNew(:) = v(:) + (h * am(1)) * g(:)
      do j = 1, nSteps
        vNew(:) = vNew(:) + (h * am(j + 1)) * history(:, j)
      end do
    end if

    ! G_{n+1} for the history, u_{n+1} as the result
    call EvaluateG(g, state, vNew, t1)
    call Propagator_Exponential_PushHistory(history, nHistory, g)

    uLast(:) = state(:)
    vLast(:) = vNew(:)
    tLast = t1
    hLast = h

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> g = G(t, v) = E(tRef - t) N(t, u) with u = E(t - tRef) v (also returned).
  subroutine EvaluateG(g, u, v, time)
    use M_Method
    use M_Propagator_Exponential

    complex(R64), intent(out), contiguous :: g(:)
    complex(R64), intent(out), contiguous :: u(:)
    complex(R64), intent(in), contiguous :: v(:)
    real(R64), intent(in) :: time

    u(:) = v(:)
    call Method_PropagateExactPart(u, tRef, time)
    call Propagator_Exponential_NonExactDerivative(g, u, time)
    call Method_PropagateExactPart(g, time, tRef)

  end subroutine

end submodule
