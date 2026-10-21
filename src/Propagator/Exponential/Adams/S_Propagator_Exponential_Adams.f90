! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Exponential Adams PECE step with the history N_n, ..., N_{n-k+1}:
!>   P: uP      = phi_0 u_n + sum_{m<k}  h phi_{m+1} xP_m,   xP = CP N
!>   E: NP      = N(t_{n+1}, uP)
!>   C: u_{n+1} = phi_0 u_n + sum_{m<=k} h phi_{m+1} xC_m,   xC = CC [NP, N]
!>   E: N_{n+1} = N(t_{n+1}, u_{n+1})
!> with phi_j = phi_j(h L). x_m are the scaled Taylor coefficients h^m p^{(m)}(t_n)
!> of the interpolation polynomial p of N (nodes θ = 0, -1, ... for the
!> predictor and θ = 1, 0, -1, ... for the corrector, in units of h); CP and
!> CC are fixed by the nodes. The history is built with Krogstad steps on the
!> first call and whenever a step does not continue the previous one. A state
!> changed by the caller between steps (absorber) is used as it is.
submodule(M_Propagator_Exponential_Adams) S_Propagator_Exponential_Adams

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Number of history values k (propagator.exponential.adams.nSteps).
  integer(I32) :: nSteps

  !> Taylor coefficients from the node values: predictor (k x k) and
  !> corrector (k + 1 x k + 1, predicted value first).
  real(R64), allocatable :: cP(:, :)
  real(R64), allocatable :: cC(:, :)

  !> History N_n, N_{n-1}, ... (newest in column 1) and its number of valid columns.
  complex(R64), allocatable :: history(:, :)
  integer(I32) :: nHistory = 0

  !> End time and step size of the last step.
  real(R64) :: tLast, hLast

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Reads `"propagator.exponential.adams.nSteps"` (1..6, default 4: order 5).
  module subroutine Propagator_Exponential_Adams_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Propagator
    use M_Propagator_Exponential

    integer(I32) :: j

    call Say_Fabricate("propagator.exponential.adams")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    nSteps = Json_Get("propagator.exponential.adams.nSteps", 4)
    if (nSteps < 1 .or. nSteps > 6) error stop "propagator.exponential.adams.nSteps must be in 1..6"

    cP = Propagator_Exponential_TaylorCoefficients([(real(1 - j, kind=R64), j=1, nSteps)])
    cC = Propagator_Exponential_TaylorCoefficients([(real(2 - j, kind=R64), j=1, nSteps + 1)])
    nHistory = 0

    Propagator_Propagate => Propagate

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Propagate(state, t0, t1)
    use M_Propagator_Exponential

    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1

    complex(R64), allocatable :: g(:), uP(:), values(:, :)
    real(R64) :: h

    h = t1 - t0
    ! (re)allocate for the first call, a new nSteps or a new state size
    if (allocated(history)) then
      if (size(history, 1) .ne. size(state) .or. size(history, 2) .ne. nSteps) deallocate (history)
    end if
    if (.not. allocated(history)) then
      allocate (history(size(state), nSteps))
      nHistory = 0
    end if
    allocate (g, uP, mold=state)

    if (nHistory > 0) then
      if (.not. Propagator_Exponential_ContinuesQ(t0, t1, tLast, hLast)) nHistory = 0
    end if

    if (nHistory .eq. 0) then
      call Propagator_Exponential_NonExactDerivative(g, state, t0)
      call Propagator_Exponential_PushHistory(history, nHistory, g)
    end if

    !------------------------------------
    ! step: Krogstad until the history is complete, then PECE
    !------------------------------------

    if (nHistory < nSteps) then
      call Propagator_Exponential_KrogstadStep(state, t0, t1)
    else
      ! predictor
      call AdamsStep(uP, state, history, cP, h)
      call Propagator_Exponential_NonExactDerivative(g, uP, t1)

      ! corrector
      allocate (values(size(state), nSteps + 1))
      values(:, 1) = g(:)
      values(:, 2:) = history(:, :)
      call AdamsStep(uP, state, values, cC, h)
      state(:) = uP(:)
    end if

    ! N_{n+1} for the history
    call Propagator_Exponential_NonExactDerivative(g, state, t1)
    call Propagator_Exponential_PushHistory(history, nHistory, g)

    tLast = t1
    hLast = h

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> uNew = phi_0(h L) u + sum_m h phi_{m+1}(h L) x_m with x = c values.
  subroutine AdamsStep(uNew, u, values, c, h)
    use M_Method

    complex(R64), intent(out), contiguous :: uNew(:)
    complex(R64), intent(in), contiguous :: u(:)
    complex(R64), intent(in), contiguous :: values(:, :)
    real(R64), intent(in) :: c(:, :)
    real(R64), intent(in) :: h

    complex(R64), allocatable :: vecs(:, :)
    integer(I32) :: m, j

    allocate (vecs(size(u), 0:size(c, 1)))
    vecs(:, 0) = u(:)
    do m = 1, size(c, 1)
      vecs(:, m) = (0.0_R64, 0.0_R64)
      do j = 1, size(c, 2)
        vecs(:, m) = vecs(:, m) + (h * c(m, j)) * values(:, j)
      end do
    end do

    call Method_ApplyExactPartPhi(uNew, vecs, h)

  end subroutine

end submodule
