! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased) S_Method_Mb_OrbBased_FieldFreeCayley
  !-----------------------------------------------------------------------------
  ! Exact field-free part h0 = T + V_static on grids without a dense h0
  ! eigen-decomposition (method.mb.orbBased.exactFieldFreeQ on non-Ylm grids
  ! whose kinetic backend binds SysKinetic_SolveShifted, e.g. grid.fem): the
  ! *Cayley procedures, which Method_Mb_OrbBased_Fabricate binds to
  ! Method_Mb_OrbBased_SetupFieldFree and _PropagateFieldFree.
  !
  ! Every stored orbital is propagated by the Crank-Nicolson (Cayley) step
  !   phi <- (1 + A)^-1 (1 - A) phi = 2 (1 + A)^-1 phi - phi,  A = i dt h0/2,
  ! which is unitary in the grid metric (h0 is self-adjoint there) for every dt
  ! and time-reversible, Cayley(-dt) = Cayley(dt)^-1. Per eigenvalue E of h0 the
  ! phase is 2 atan(E dt/2) instead of E dt: second order (local error O(dt^3)),
  ! like the Strang composition it is used in (propagator.splitStep with
  ! exactPartScheme strang).
  !
  ! The system (1 + A) y = phi is solved by the flexible GMRES(m) of
  ! M_Utils_Gmres in the grid metric (Grid_InnerProduct), right-preconditioned
  ! with the exact inverse of its kinetic factor,
  !   P^-1 = (1 + i dt T/2)^-1 = sigma (T + sigma)^-1,  sigma = 2/(i dt),
  ! (SysKinetic_SolveShifted with a complex shift; fast diagonalization on
  ! grid.fem), so that the Krylov iteration only sees the potential,
  ! (1 + A) P^-1 = 1 + (i dt/2) V P^-1: a few iterations per step also with the
  ! Galerkin nuclear spikes. This file supplies the physics as the three
  ! callbacks of the solver: ApplyPreconditioner (P^-1), the Arnoldi operator
  ! ApplyPreconditionedSystem, v + (i dt/2) V z for z = P^-1 v (one
  ! preconditioner solve and one application of V per iteration, no kinetic
  ! operator), and ApplySystem, (1 + A) y with h0 applied explicitly
  ! (SysKinetic_MultiplyWithKineticOp + V_static) for the true residuals.
  ! V_static is always applied through the SysPotential contract
  ! (SysPotential_MultiplyWithExternalPotential; the backend owns its
  ! representation, a nodal product on grid.fem). Initial guess y0 = P^-1 phi;
  ! every convergence decision, the first one included, is made on the true
  ! residual b - (1 + A) y, which absorbs the round-off of P P^-1 = 1 (for y0
  ! it equals -(i dt/2) V y0, with the mean-field reference below
  ! -(i dt/2) (V + U) y0, up to that round-off): ||b - (1 + A) y|| <=
  ! exactFieldFreeTol ||b|| (grid metric), at most 1 + exactFieldFreeMaxRestarts
  ! cycles of exactFieldFreeKrylovDim iterations, error stop otherwise.
  !
  ! The callbacks read the step data stepDt and stepBt (set per propagation
  ! and per orbital) and the scratch vector. The GMRES context (2m + 2 work
  ! vectors) and the solution and scratch vectors are module data, allocated
  ! at the first propagation after a setup and kept (their pages stay mapped
  ! from call to call), so Method_Mb_OrbBased_PropagateFieldFreeCayley is not
  ! reentrant: it stops inside an active OpenMP parallel region (the methods
  ! call it serially; its vector passes and the kinetic solves are parallel
  ! inside).
  !
  ! Not available (checked at fabricate or here): phi functions (MCTDHX binds
  ! no Method_ApplyExactPartPhi, so propagator.exponential stops), the spectral
  ! cut, imaginary time (Cayley does not damp the high modes).
  !
  ! Frozen mean-field reference (method.mb.orbBased.exactFieldFreeMeanFieldQ,
  ! MCTDHX; S_Method_Mb_OrbBased_FieldFreeMeanField): once it is built, h0 is
  ! replaced by A = h0 + U with a static Hermitian one-body operator U. The
  ! Arnoldi operator and the true residual apply A, U per vector through
  ! Method_Mb_OrbBased_ApplyFieldFreeMeanFieldOnOrb (a no-op without the
  ! reference); the preconditioner is unchanged. The setup ends with
  ! Method_Mb_OrbBased_SetupFieldFreeMeanField, which releases U.
  !-----------------------------------------------------------------------------

  use M_Utils_Gmres

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Whether Method_Mb_OrbBased_SetupFieldFreeCayley has run.
  logical :: setupQ = .false.

  !> Whether a static external potential exists ("sysPotential"); without it
  !> V_static = 0 and the SysPotential procedures are never called.
  logical :: potentialQ = .false.

  !> Static external potential V_static = SysPotential_FillExternalPotential
  !> at time 0 in the representation of the potential backend, cached in the
  !> setup (allocated only with potentialQ) and applied with
  !> SysPotential_MultiplyWithExternalPotential.
  complex(R64), allocatable :: staticPotential(:)

  !> GMRES context (M_Utils_Gmres, Gmres_CreateCtx in the setup): the controls
  !> method.mb.orbBased.exactFieldFreeTol (relative residual in the grid
  !> metric), ...KrylovDim (restart length m) and ...MaxRestarts, and the
  !> solver workspace (2m + 2 vectors of nPoints, allocated by the first solve
  !> after the setup and kept).
  type(T_Gmres_Ctx) :: gmresCtx

  !> Step data of the solver callbacks: the time step of the propagation and
  !> the body type of the orbital being solved for.
  real(R64) :: stepDt = 0.0_R64
  integer(I32) :: stepBt = 1

  !> Work vectors of the propagation (AllocateWorkArrays): the solution and
  !> a scratch vector of the callbacks (U or T times a vector).
  complex(R64), allocatable :: workSolution(:)
  complex(R64), allocatable :: workScratch(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_SetupFieldFreeCayley
    !---------------------------------------------------------------------------
    ! Requirements of a static h0 (time-independent kinetic operator, time-
    ! and body-type-independent potential), V_static and the GMRES controls.
    ! The kinetic operator may depend on the body type (masses): every orbital
    ! uses the operators of its body type. Creates the GMRES context and
    ! releases the work vectors of an earlier setup (the next propagation
    ! allocates them for this one) and ends with the setup of the frozen
    ! mean-field reference (Method_Mb_OrbBased_SetupFieldFreeMeanField:
    ! releases it, checks the interaction, prints its line after the GMRES
    ! line).
    !---------------------------------------------------------------------------
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysPotential

    real(R64) :: tolerance
    integer(I32) :: krylovDim, maxRestarts

    call Say_Setup("method.mb.orbBased.exactFieldFree (Cayley)")

    !------------------------------------
    ! applicability
    !------------------------------------

    if (.not. associated(SysKinetic_SolveShifted)) &
      error stop "method.mb.orbBased.exactFieldFreeQ (Cayley): the kinetic backend does not bind SysKinetic_SolveShifted"
    if (.not. SysKinetic_timeIndependentQ) &
      error stop "method.mb.orbBased.exactFieldFreeQ requires a time-independent kinetic operator"
    if (Json_Get("method.mb.orbBased.exactFieldFreeGpuQ", .false.)) &
      error stop "method.mb.orbBased.exactFieldFreeGpuQ is only available on Ylm grids (dense h0 blocks)"

    !------------------------------------
    ! static potential
    !------------------------------------

    if (allocated(staticPotential)) deallocate (staticPotential)
    potentialQ = Json_GetExistence("sysPotential")
    if (potentialQ) then
      if (.not. (SysPotential_timeIndependentQ .and. SysPotential_bodyTypeIndependentQ)) &
        error stop "method.mb.orbBased.exactFieldFreeQ requires a time- and body-type-independent potential &
                   &(a laser belongs to sysGauge, e.g. sysGauge.fem.length)"
      call SysPotential_FillExternalPotential(staticPotential, 0.0_R64)
    end if

    !------------------------------------
    ! GMRES controls
    !------------------------------------

    tolerance = Json_Get("method.mb.orbBased.exactFieldFreeTol", 1.0e-12_R64)
    krylovDim = Json_Get("method.mb.orbBased.exactFieldFreeKrylovDim", 30)
    maxRestarts = Json_Get("method.mb.orbBased.exactFieldFreeMaxRestarts", 20)

    if (.not. (tolerance > 0.0_R64)) error stop "method.mb.orbBased.exactFieldFreeTol must be > 0"
    if (krylovDim < 1) error stop "method.mb.orbBased.exactFieldFreeKrylovDim must be >= 1"
    if (maxRestarts < 0) error stop "method.mb.orbBased.exactFieldFreeMaxRestarts must be >= 0"

    call Gmres_CreateCtx(gmresCtx, krylovDim, maxRestarts, tolerance)
    call ReleaseWorkArrays
    setupQ = .true.

    Method_Mb_OrbBased_cayleyNSolves = 0
    Method_Mb_OrbBased_cayleyNIterations = 0
    Method_Mb_OrbBased_cayleyMaxIterations = 0

    write (*, '(A,ES9.2,A,I0,A,I0)') " exact field-free part: Cayley step, GMRES tolerance", tolerance, &
      ", Krylov dimension ", krylovDim, ", restarts ", maxRestarts

    !------------------------------------
    ! frozen mean-field reference
    !------------------------------------

    call Method_Mb_OrbBased_SetupFieldFreeMeanField

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_PropagateFieldFreeCayley(orbs, dt, imagTimeQ)
    !---------------------------------------------------------------------------
    ! phi <- 2 y - phi with (1 + i dt A/2) y = phi, A = h0 (+ U with the
    ! frozen mean-field reference), per stored orbital (body type of the
    ! orbital; restricted storage: body type 1): Gmres_Solve with the
    ! callbacks below, which read the step data stepDt and stepBt; not
    ! converged: message and error stop. The orbital loop is serial and the
    ! context, the work vectors and the step data are shared (not reentrant),
    ! the vector passes and the kinetic solves are parallel internally.
    !---------------------------------------------------------------------------
    !$ use omp_lib, only: omp_in_parallel
    use M_Grid

    complex(R64), intent(inout), contiguous :: orbs(:, :)
    real(R64), intent(in) :: dt
    logical, intent(in) :: imagTimeQ

    real(R64) :: relativeResidual
    integer(I32) :: iOrb, nIterations
    logical :: inParallelQ, convergedQ

    if (imagTimeQ) &
      error stop "method.mb.orbBased.exactFieldFreeQ: the Cayley exact part (non-Ylm grids) is real-time only &
                 &(it does not damp); use the SCF ground state (groundSolver) instead of imaginary time"
    if (.not. setupQ) &
      error stop "method.mb.orbBased.exactFieldFreeQ (Cayley): call Method_Setup before propagating"
    inParallelQ = .false.
    !$ inParallelQ = omp_in_parallel()
    if (inParallelQ) &
      error stop "method.mb.orbBased.exactFieldFreeQ (Cayley): the propagation is not reentrant &
                 &(shared work arrays); call it outside OpenMP parallel regions"
    if (abs(dt) <= 0.0_R64) return

    call AllocateWorkArrays(Grid_nPoints)
    stepDt = dt

    do iOrb = 1, size(orbs, 2)
      stepBt = Method_Mb_OrbBased_bodyTypeOfOrb(iOrb)
      call Gmres_Solve(workSolution, nIterations, relativeResidual, convergedQ, gmresCtx, orbs(:, iOrb), &
                       ApplyPreconditioner, ApplyPreconditionedSystem, ApplySystem, Grid_InnerProduct)
      if (.not. convergedQ) then
        write (*, '(A,ES10.3,A,ES9.2,A,I0,A,I0,A)') " Cayley GMRES: relative residual", relativeResidual, &
          " above the tolerance", gmresCtx % tolerance, " after ", nIterations, " iterations (", &
          gmresCtx % maxRestarts, " restarts)"
        error stop "method.mb.orbBased.exactFieldFreeQ: GMRES of the Cayley step did not converge &
                   &(raise exactFieldFreeKrylovDim/exactFieldFreeMaxRestarts or exactFieldFreeTol, or lower dt)"
      end if
      call CountSolve(nIterations)
      call CombineCayley(orbs(:, iOrb), workSolution)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Allocates the work vectors for nG points unless they exist with this
  !> size (the GMRES workspace is allocated by its solve).
  subroutine AllocateWorkArrays(nG)
    integer(I32), intent(in) :: nG

    if (allocated(workSolution)) then
      if (size(workSolution) .eq. nG) return
      call ReleaseWorkArrays
    end if

    allocate (workSolution(nG), workScratch(nG))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Deallocates the work vectors.
  subroutine ReleaseWorkArrays

    if (allocated(workSolution)) deallocate (workSolution)
    if (allocated(workScratch)) deallocate (workScratch)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> z = P^-1 v = sigma (T + sigma)^-1 v, sigma = 2/(i dt) = -2i/dt (the
  !> inverse of the kinetic factor 1 + i dt T/2 of body type bt): the
  !> preconditioner of the GMRES solve (I_Gmres_ApplyPrecOnVec; dt = stepDt,
  !> bt = stepBt).
  subroutine ApplyPreconditioner(z, v)
    !$ use omp_lib, only: omp_in_parallel
    use M_SysKinetic

    complex(R64), intent(out), contiguous :: z(:)
    complex(R64), intent(in), contiguous :: v(:)

    complex(R64) :: sigma
    integer(I32) :: i

    sigma = cmplx(0.0_R64, -2.0_R64 / stepDt, kind=R64)
    call SysKinetic_SolveShifted(z, v, sigma, bt_=stepBt)

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(z)
      z(i) = sigma * z(i)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> w = V_static z through the SysPotential contract
  !> (SysPotential_MultiplyWithExternalPotential: the backend applies its own
  !> representation of the cached potential, a nodal product on grid.fem);
  !> w = 0 without sysPotential.
  subroutine ApplyStaticPotential(w, z)
    !$ use omp_lib, only: omp_in_parallel
    use M_SysPotential

    complex(R64), intent(out), contiguous :: w(:)
    complex(R64), intent(in), contiguous :: z(:)

    integer(I32) :: i

    if (potentialQ) then
      call SysPotential_MultiplyWithExternalPotential(w, staticPotential, z)
      return
    end if

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(w)
      w(i) = (0.0_R64, 0.0_R64)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> w = (V_static + U) z: the static potential (ApplyStaticPotential) plus,
  !> once the frozen mean-field reference is built and active, U z of body
  !> type bt (Method_Mb_OrbBased_ApplyFieldFreeMeanFieldOnOrb; scratch: its
  !> work vector).
  subroutine ApplyStaticOperators(w, scratch, z, bt)
    complex(R64), intent(out), contiguous :: w(:)
    complex(R64), intent(inout), contiguous :: scratch(:)
    complex(R64), intent(in), contiguous :: z(:)
    integer(I32), intent(in) :: bt

    call ApplyStaticPotential(w, z)
    call Method_Mb_OrbBased_ApplyFieldFreeMeanFieldOnOrb(w, scratch, z, bt)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> w = (1 + A) P^-1 v = v + (i dt/2) (V_static + U) z for z = P^-1 v: the
  !> Arnoldi operator of the GMRES solve (I_Gmres_ApplyMatOnPrecVec; dt =
  !> stepDt, body type stepBt, workScratch: work vector of U).
  subroutine ApplyPreconditionedSystem(w, v, z)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(out), contiguous :: w(:)
    complex(R64), intent(in), contiguous :: v(:)
    complex(R64), intent(in), contiguous :: z(:)

    complex(R64) :: halfStep
    integer(I32) :: i

    halfStep = cmplx(0.0_R64, 0.5_R64 * stepDt, kind=R64)
    call ApplyStaticOperators(w, workScratch, z, stepBt)

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(w)
      w(i) = v(i) + halfStep * w(i)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> w = (1 + i dt A/2) z with A z = T z (body type bt) + V_static z (+ U z):
  !> the explicit operator of the true residuals of the GMRES solve
  !> (I_Gmres_ApplyMatOnVec; dt = stepDt, bt = stepBt, workScratch: work
  !> vector of U, then T z).
  subroutine ApplySystem(w, z)
    !$ use omp_lib, only: omp_in_parallel
    use M_SysKinetic

    complex(R64), intent(out), contiguous :: w(:)
    complex(R64), intent(in), contiguous :: z(:)

    complex(R64) :: halfStep
    integer(I32) :: i

    halfStep = cmplx(0.0_R64, 0.5_R64 * stepDt, kind=R64)
    call ApplyStaticOperators(w, workScratch, z, stepBt)
    call SysKinetic_MultiplyWithKineticOp(workScratch, z, 0.0_R64, bt_=stepBt)

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(w)
      w(i) = z(i) + halfStep * (workScratch(i) + w(i))
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> phi <- 2 y - phi (Cayley transform from the solution of (1 + A) y = phi).
  subroutine CombineCayley(phi, y)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(inout), contiguous :: phi(:)
    complex(R64), intent(in), contiguous :: y(:)

    integer(I32) :: i

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(phi)
      phi(i) = 2.0_R64 * y(i) - phi(i)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Records one GMRES solve with nIterations iterations.
  subroutine CountSolve(nIterations)
    integer(I32), intent(in) :: nIterations

    Method_Mb_OrbBased_cayleyNSolves = Method_Mb_OrbBased_cayleyNSolves + 1
    Method_Mb_OrbBased_cayleyNIterations = Method_Mb_OrbBased_cayleyNIterations + nIterations
    Method_Mb_OrbBased_cayleyMaxIterations = max(Method_Mb_OrbBased_cayleyMaxIterations, nIterations)

  end subroutine

end submodule
