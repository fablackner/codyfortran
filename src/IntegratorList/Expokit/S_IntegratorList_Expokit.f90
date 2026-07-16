! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_IntegratorList_Expokit.f90
!> @brief Krylov subspace matrix exponential integrator (Expokit-style).
!>
!> @details
!> Computes the action of the matrix exponential on a vector:
!>
!>     ψ(t + Δt) = exp(−i Δt Ĥ) ψ(t)
!>
!> using an Arnoldi/Lanczos-based Krylov subspace approximation. This approach
!> is highly efficient for large, sparse Hamiltonians where forming the full
!> exponential is infeasible.
!>
!> ## Algorithm Outline
!>
!> 1. Build a Krylov basis V_m spanning {v, Av, A²v, …, A^{m-1}v}.
!> 2. Project A onto V_m to get a small Hessenberg matrix H_m.
!> 3. Compute exp(Δt · H_m) via dense methods.
!> 4. Map back to the full space: ψ_new ≈ ||v|| · V_m · exp(Δt · H_m) · e₁.
!>
!> ## Properties
!>
!> - **Accuracy**: Controlled by `krylovDim` and `tolerance`.
!> - **Stability**: Exact for linear time-independent Hamiltonians.
!> - **Cost**: O(m) matrix-vector products per step; m ≪ N.
!>
!> ## Configuration (JSON)
!>
!> | Key         | Type    | Default | Description                     |
!> |-------------|---------|---------|---------------------------------|
!> | `krylovDim` | integer | 30      | Krylov subspace dimension       |
!> | `tolerance` | real    | 1e-7    | Error tolerance for exp approx  |
!> | `maxSteps`  | integer | 1000    | Max internal sub-steps          |
!>
!> @note The `TimeDerivative` callback supplies −i Ĥ |ψ⟩. Internally we
!>   multiply by i to recover Ĥ |ψ⟩ for the Expokit library.
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

    this % krylovDim = Json_Get("krylovDim", 30, path_=this % path)
    this % tolerance = Json_Get("tolerance", 1.0e-7_R64, path_=this % path)
    this % maxSteps = Json_Get("maxSteps", 1000, path_=this % path)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Expokit), intent(inout) :: this

    call Say_Setup(this % path//".expokit")

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Advance the state via Krylov-based matrix exponential.
  !>
  !> Computes ψ(t1) ≈ exp(−i (t1−t0) Ĥ) ψ(t0) using the Expokit library's
  !> `zhexpv` routine for Hermitian matrices.
  !>
  !> @param[inout] this   The Expokit integrator instance.
  !> @param[inout] state  Complex state vector (modified in-place).
  !> @param[in]    t0     Starting time.
  !> @param[in]    t1     Target ending time.
  module subroutine Integrate(this, state, t0, t1)
    use M_Utils_ExpokitLib, only: ExpokitLib_IntegrateSym
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Expokit), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1

    integer :: iFlag, nSteps
    real(R64) :: tStep
    character(len=100) :: errorMsg

    tStep = t1 - t0
    nSteps = this % maxSteps

    ! Call Expokit wrapper for Hermitian/symmetric matrices
    call ExpokitLib_IntegrateSym(state, tStep, ApplyMatrix, &
                                 this % krylovDim, this % tolerance, nSteps, iFlag)

    if (iFlag .ne. 0) then
      write (errorMsg, '(A,I0)') "Expokit zhexpv failed with status: ", iFlag
      error stop errorMsg
    end if

  contains
    !> Internal callback: convert TimeDerivative output to H*x for Expokit.
    !>
    !> TimeDerivative computes −i Ĥ x; Expokit needs Ĥ x, so we multiply by i.
    subroutine ApplyMatrix(dState, inState)
      complex(R64), intent(out), contiguous, target :: dState(:)
      complex(R64), intent(in), contiguous, target  :: inState(:)

      call this % TimeDerivative(dState, inState, t0)
      dState = (0.0_R64, 1.0_R64) * dState  ! Multiply by i: −i·(−i Ĥ x) = Ĥ x
    end subroutine
  end subroutine

end submodule
