! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Factory of the exponential integrators and their shared building blocks.
!>
!> Krogstad step, phi_j = phi_j(h L), phi_j' = phi_j(h L / 2):
!>   a       = phi_0' u + h/2 phi_1' N_u
!>   b       = phi_0' u + h/2 phi_1' N_u + h phi_2' (N_a - N_u)
!>   c       = phi_0 u + h phi_1 N_u + 2 h phi_2 (N_b - N_u)
!>   u_{n+1} = phi_0 u + h phi_1 N_u + h phi_2 (-3 N_u + 2 N_a + 2 N_b - N_c)
!>                     + 4 h phi_3 (N_u - N_a - N_b + N_c)
!> with N_x = N(t_x, x) at t_a = t_b = t_n + h/2 and t_c = t_n + h.
submodule(M_Propagator_Exponential) S_Propagator_Exponential

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Exactly one scheme key must be present:
!>   - `"propagator.exponential.krogstad4"` → Krogstad step (this submodule)
!>   - `"propagator.exponential.adams"` → Propagator_Exponential_Adams_Fabricate
!>   - `"propagator.exponential.adamsLawson"` → Propagator_Exponential_AdamsLawson_Fabricate
  module subroutine Propagator_Exponential_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_IntegratorList
    use M_Method
    use M_Propagator
    use M_Propagator_Exponential_Adams
    use M_Propagator_Exponential_AdamsLawson

    call Say_Fabricate("propagator.exponential")

    !------------------------------------
    ! applicability
    !------------------------------------

    ! The Cayley exact part (method.mb.orbBased.exactFieldFreeQ on non-Ylm
    ! grids, e.g. grid.fem) binds Method_PropagateExactPart but no phi functions
    if (.not. (associated(Method_PropagateExactPart) .and. associated(Method_ApplyExactPartPhi))) &
      error stop "propagator.exponential requires a method with an exactly propagated part and its phi functions &
                 &(method.mb.orbBased.exactFieldFreeQ with mctdhx on a Ylm grid); exactFieldFreeQ on non-Ylm grids &
                 &(e.g. grid.fem) is a Cayley step without phi functions: use propagator.splitStep with &
                 &exactPartScheme strang"
    if (IntegratorList(1) % e % imagTimeQ .or. IntegratorList(2) % e % imagTimeQ) &
      error stop "propagator.exponential: imaginary time is not supported"

    !------------------------------------
    ! branch: dispatch to the scheme
    !------------------------------------

    if (Json_GetExistence("propagator.exponential.krogstad4")) then
      Propagator_Propagate => Propagator_Exponential_KrogstadStep

    else if (Json_GetExistence("propagator.exponential.adams")) then
      call Propagator_Exponential_Adams_Fabricate

    else if (Json_GetExistence("propagator.exponential.adamsLawson")) then
      call Propagator_Exponential_AdamsLawson_Fabricate

    else
      error stop "propagator.exponential is missing one of: krogstad4, adams, adamsLawson"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Propagator_Exponential_KrogstadStep(state, t0, t1)
    use M_Method

    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1

    complex(R64), allocatable :: nU(:), nA(:), nB(:), nC(:), stage(:), vecs(:, :)
    real(R64) :: h, tHalf

    h = t1 - t0
    tHalf = t0 + 0.5_R64 * h
    allocate (nU, nA, nB, nC, stage, mold=state)
    allocate (vecs(size(state), 4))

    call Propagator_Exponential_NonExactDerivative(nU, state, t0)

    vecs(:, 1) = state(:)
    vecs(:, 2) = (0.5_R64 * h) * nU(:)
    call Method_ApplyExactPartPhi(stage, vecs(:, 1:2), 0.5_R64 * h)
    call Propagator_Exponential_NonExactDerivative(nA, stage, tHalf)

    vecs(:, 3) = h * (nA(:) - nU(:))
    call Method_ApplyExactPartPhi(stage, vecs(:, 1:3), 0.5_R64 * h)
    call Propagator_Exponential_NonExactDerivative(nB, stage, tHalf)

    vecs(:, 2) = h * nU(:)
    vecs(:, 3) = (2.0_R64 * h) * (nB(:) - nU(:))
    call Method_ApplyExactPartPhi(stage, vecs(:, 1:3), h)
    call Propagator_Exponential_NonExactDerivative(nC, stage, t1)

    vecs(:, 3) = h * (-3.0_R64 * nU(:) + 2.0_R64 * nA(:) + 2.0_R64 * nB(:) - nC(:))
    vecs(:, 4) = (4.0_R64 * h) * (nU(:) - nA(:) - nB(:) + nC(:))
    call Method_ApplyExactPartPhi(state, vecs, h)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Propagator_Exponential_NonExactDerivative(dState, state, time)
    use M_IntegratorList

    complex(R64), intent(out), contiguous :: dState(:)
    complex(R64), intent(in), contiguous :: state(:)
    real(R64), intent(in) :: time

    complex(R64), allocatable :: d2(:)

    allocate (d2, mold=state)
    call IntegratorList(1) % e % TimeDerivative(dState, state, time)
    call IntegratorList(2) % e % TimeDerivative(d2, state, time)
    dState(:) = dState(:) + d2(:)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> C = diag(m!) V^{-1} with the Vandermonde matrix V_{jm} = θ_j^m.
  module function Propagator_Exponential_TaylorCoefficients(nodes) result(res)
    use M_Utils_LapackLib

    real(R64), intent(in) :: nodes(:)
    real(R64), allocatable :: res(:, :)

    real(R64), allocatable :: vandermonde(:, :)
    integer(I32) :: n, j, m

    n = size(nodes)
    allocate (vandermonde(n, n), res(n, n))
    res(:, :) = 0.0_R64
    do j = 1, n
      res(j, j) = 1.0_R64
      do m = 0, n - 1
        vandermonde(j, m + 1) = nodes(j)**m
      end do
    end do

    call LapackLib_SolveLinearSystem(vandermonde, res)

    do m = 1, n - 1
      res(m + 1, :) = res(m + 1, :) * gamma(real(m + 1, kind=R64))
    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Propagator_Exponential_PushHistory(history, nHistory, g)
    complex(R64), intent(inout), contiguous :: history(:, :)
    integer(I32), intent(inout) :: nHistory
    complex(R64), intent(in), contiguous :: g(:)

    integer(I32) :: j

    nHistory = min(nHistory + 1, size(history, 2))
    do j = nHistory, 2, -1
      history(:, j) = history(:, j - 1)
    end do
    history(:, 1) = g(:)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function Propagator_Exponential_ContinuesQ(t0, t1, tLast, hLast) result(res)
    real(R64), intent(in) :: t0, t1
    real(R64), intent(in) :: tLast, hLast
    logical :: res

    real(R64) :: tol

    tol = 1e-8_R64 * abs(t1 - t0)
    res = abs(t0 - tLast) <= tol .and. abs((t1 - t0) - hLast) <= tol

  end function

end submodule
