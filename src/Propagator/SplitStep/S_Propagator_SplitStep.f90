! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Factory for split-operator propagator backends.
!>
!> Dispatches to the appropriate composition order:
!> - `order2`: Symmetric Strang splitting, O(Δt³) local error
!> - `order4`: Yoshida/Forest–Ruth composition, O(Δt⁵) local error
!>
!> The split-step method decomposes Ĥ = Â + B̂ and approximates:
!>   exp(−iĤΔt) ≈ product of exp(−iÂ·aₖΔt) and exp(−iB̂·bₖΔt)
!>
!> Each operator (Â, B̂) is applied via a separate integrator in IntegratorList.
!> Every Â step goes through Propagator_SplitStep_IntegrateLin, which composes
!> it with Method_PropagateExactPart when the method binds one
!> (propagator.splitStep.exactPartScheme), and otherwise calls IntegratorList(1)
!> alone.
submodule(M_Propagator_SplitStep) S_Propagator_SplitStep

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Composition of the linear step when the method propagates part of it
  !> exactly (`Method_PropagateExactPart`): exact substeps exactCoeffs(k) Δt
  !> alternate with integrator-1 substeps integratorCoeffs(k) Δt,
  !> exact(1) int(1) exact(2) ... int(n) exact(n+1)
  !> (propagator.splitStep.exactPartScheme).
  real(R64), allocatable :: exactCoeffs(:)
  real(R64), allocatable :: integratorCoeffs(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Reads `"propagator.splitStep"` and delegates to order-specific backend.
!>
!> Exactly one order key must be present:
!>   - `"propagator.splitStep.order2"` → 2nd-order Strang splitting
!>   - `"propagator.splitStep.order4"` → 4th-order Yoshida composition
!>
!> Optional: `"propagator.splitStep.exactPartScheme"` (strang | yoshida4 |
!> blanesMoan6, default strang). It is validated always but used only if
!> Method_PropagateExactPart is associated.
  module subroutine Propagator_SplitStep_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Propagator_SplitStep_Order2
    use M_Propagator_SplitStep_Order4

    character(len=:), allocatable :: exactPartScheme
    real(R64) :: w0, w1, a1, a2, a3, b1, b2

    call Say_Fabricate("propagator.splitStep")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Splitting of the exactly propagated part and integrator 1 inside the
    ! linear step (only used if the method binds Method_PropagateExactPart)
    exactPartScheme = Json_Get("propagator.splitStep.exactPartScheme", "strang")

    select case (exactPartScheme)
    case ("strang")
      ! second order, O(Δt³) local error
      exactCoeffs = [0.5_R64, 0.5_R64]
      integratorCoeffs = [1.0_R64]

    case ("yoshida4")
      ! fourth order (Yoshida 1990), 4 exact and 3 integrator substeps
      w1 = 1.0_R64 / (2.0_R64 - 2.0_R64**(1.0_R64 / 3.0_R64))
      w0 = 1.0_R64 - 2.0_R64 * w1
      exactCoeffs = [0.5_R64 * w1, 0.5_R64 * (w0 + w1), 0.5_R64 * (w0 + w1), 0.5_R64 * w1]
      integratorCoeffs = [w1, w0, w1]

    case ("blanesMoan6")
      ! fourth order with small error constants (Blanes & Moan 2002, PRK S6),
      ! 7 exact and 6 integrator substeps
      a1 = 0.0792036964311957_R64
      a2 = 0.353172906049774_R64
      a3 = -0.0420650803577195_R64
      b1 = 0.209515106613362_R64
      b2 = -0.143851773179818_R64
      exactCoeffs = [a1, a2, a3, 1.0_R64 - 2.0_R64 * (a1 + a2 + a3), a3, a2, a1]
      integratorCoeffs = [b1, b2, 0.5_R64 - b1 - b2, 0.5_R64 - b1 - b2, b2, b1]

    case default
      error stop "propagator.splitStep.exactPartScheme must be one of: strang, yoshida4, blanesMoan6"
    end select

    !------------------------------------
    ! branch: dispatch to order-specific backend
    !------------------------------------

    if (Json_GetExistence("propagator.splitStep.order2")) then
      call Propagator_SplitStep_Order2_Fabricate

    else if (Json_GetExistence("propagator.splitStep.order4")) then
      call Propagator_SplitStep_Order4_Fabricate

    else
      error stop "propagator.splitStep is missing one of: order2, order4"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Propagator_SplitStep_IntegrateLin(state, t0, t1)
    use M_IntegratorList
    use M_Method

    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1

    real(R64) :: dt, tExact, tHalf
    integer(I32) :: k

    if (.not. associated(Method_PropagateExactPart)) then
      call IntegratorList(1) % e % Integrate(state, t0, t1)
      return
    end if

    dt = t1 - t0

    ! The exact substeps carry the time; each integrator substep is centred on
    ! the time reached so that its midpoint evaluation of the time-dependent
    ! part is consistent with the symmetric composition
    tExact = t0
    do k = 1, size(integratorCoeffs)
      call Method_PropagateExactPart(state, tExact, tExact + exactCoeffs(k) * dt)
      tExact = tExact + exactCoeffs(k) * dt

      tHalf = 0.5_R64 * integratorCoeffs(k) * dt
      call IntegratorList(1) % e % Integrate(state, tExact - tHalf, tExact + tHalf)
    end do
    call Method_PropagateExactPart(state, tExact, t1)

  end subroutine

end submodule
