!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The Cayley exact part of grid.fem (method.mb.orbBased.exactFieldFreeQ on a
!> non-Ylm grid, S_Method_Mb_OrbBased_FieldFreeCayley) on its own, including
!> the restart path of its GMRES. H atom on the mesh of
!> T_Method_Sb_LengthGaugeFem (box ±6, nLocals 4, hMin 0.5, N = 12167),
!> h0 = T + V_nuc, two complex test orbitals that are not eigenstates of h0
!> (Slater functions with plane-wave phases, normalized in the grid metric).
!>
!> For every time step of program.timeSteps (the substep dt/4 = 0.0025 of the
!> Fem dynamics tests and the contract step 0.05 of T_Fem_05) the step is
!> computed twice: with exactFieldFreeKrylovDim 30 (one GMRES cycle) and,
!> after a new setup, with exactFieldFreeKrylovDim = program.restartKrylovDim
!> = 2 and exactFieldFreeMaxRestarts = program.restartMaxRestarts = 50 (many
!> cycles, each restarted from the true residual).
!>
!> Checks of the restarted step:
!> - Crank-Nicolson residual |(1 + A) phi_new - (1 - A) phi| / |phi| with
!>   A = i dt h0/2 and h0 applied here (1e-11; GMRES tolerance 1e-12)
!> - unitarity in the grid metric: the overlap matrix of the two orbitals is
!>   unchanged (1e-11)
!> - Cayley(-dt) Cayley(dt) = 1 (relative 1e-11)
!> - agreement with the one-cycle step (relative 1e-11)
!> - the restart path ran: more GMRES iterations per solve than the Krylov
!>   dimension
!> Both runs: GMRES iterations per solve at most program.maxIterationsOneCycle
!> and program.maxIterationsRestart (measured + 2). The true-residual restarts
!> absorb a wrong Krylov operator, a wrong initial residual or a restart from
!> a stale residual norm only with many more iterations, or not at all.
program T_Method_FieldFreeCayleyFem
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Fem
  use M_SysKinetic
  use M_SysPotential
  use M_Method
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  implicit none

  !> Test orbitals
  integer(I32), parameter :: N_ORBS = 2

  type(error_type), allocatable :: error
  complex(R64), allocatable :: phi(:, :), phiOneCycle(:, :, :), phiRestart(:, :), phiBack(:, :), potential(:)
  real(R64), allocatable :: timeSteps(:), residual(:), unitarityDev(:), roundTrip(:), agreement(:)
  integer(I32), allocatable :: iterationsOneCycle(:), iterationsRestart(:), boundOneCycle(:), boundRestart(:)
  integer(I32) :: restartKrylovDim, restartMaxRestarts, nDt, iDt, iOrb
  logical :: foundQ

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Method/T_Method_FieldFreeCayleyFem.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call Method_Fabricate

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup

  call Say_Fabricate("program")
  timeSteps = Json_Get("program.timeSteps", [0.0025_R64, 0.05_R64])
  restartKrylovDim = Json_Get("program.restartKrylovDim", 2)
  restartMaxRestarts = Json_Get("program.restartMaxRestarts", 50)
  nDt = size(timeSteps)
  boundOneCycle = Json_Get("program.maxIterationsOneCycle", [(0, iDt=1, nDt)])
  boundRestart = Json_Get("program.maxIterationsRestart", [(0, iDt=1, nDt)])
  if (size(boundOneCycle) .ne. nDt .or. size(boundRestart) .ne. nDt) &
    error stop "T_Method_FieldFreeCayleyFem: one iteration bound per time step needed"

  !==========================================
  call Say_Section("perform action")
  !==========================================

  call FillTestOrbitals(phi)
  call SysPotential_FillExternalPotential(potential, 0.0_R64)
  allocate (phiOneCycle(Grid_nPoints, N_ORBS, nDt), phiRestart(Grid_nPoints, N_ORBS), phiBack(Grid_nPoints, N_ORBS))
  allocate (iterationsOneCycle(nDt), iterationsRestart(nDt))
  allocate (residual(nDt), unitarityDev(nDt), roundTrip(nDt), agreement(nDt))

  ! One cycle: Krylov dimension 30 of the JSON
  call Method_Mb_OrbBased_SetupFieldFree
  do iDt = 1, nDt
    phiOneCycle(:, :, iDt) = phi(:, :)
    call PropagateCounted(iterationsOneCycle(iDt), phiOneCycle(:, :, iDt), timeSteps(iDt))
  end do

  ! Restarted: Krylov dimension 2, new setup
  call jsonFile % update("method.mb.orbBased.exactFieldFreeKrylovDim", restartKrylovDim, foundQ)
  if (.not. foundQ) error stop "T_Method_FieldFreeCayleyFem: exactFieldFreeKrylovDim missing in the JSON"
  call jsonFile % update("method.mb.orbBased.exactFieldFreeMaxRestarts", restartMaxRestarts, foundQ)
  if (.not. foundQ) error stop "T_Method_FieldFreeCayleyFem: exactFieldFreeMaxRestarts missing in the JSON"
  call Method_Mb_OrbBased_SetupFieldFree

  do iDt = 1, nDt
    phiRestart(:, :) = phi(:, :)
    call PropagateCounted(iterationsRestart(iDt), phiRestart, timeSteps(iDt))

    residual(iDt) = 0.0_R64
    agreement(iDt) = 0.0_R64
    do iOrb = 1, N_ORBS
      residual(iDt) = max(residual(iDt), CrankNicolsonResidual(phiRestart(:, iOrb), phi(:, iOrb), timeSteps(iDt)))
      agreement(iDt) = max(agreement(iDt), MetricNorm(phiRestart(:, iOrb) - phiOneCycle(:, iOrb, iDt)))
    end do
    unitarityDev(iDt) = maxval(abs(OverlapMatrix(phiRestart) - OverlapMatrix(phi)))

    phiBack(:, :) = phiRestart(:, :)
    call Method_Mb_OrbBased_PropagateFieldFree(phiBack, -timeSteps(iDt), .false.)
    roundTrip(iDt) = 0.0_R64
    do iOrb = 1, N_ORBS
      roundTrip(iDt) = max(roundTrip(iDt), MetricNorm(phiBack(:, iOrb) - phi(:, iOrb)))
    end do
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(1X, A, I0, A, I0, A, I0)') "mesh N = ", Grid_nPoints, ", restarted runs: Krylov dimension ", &
    restartKrylovDim, ", restarts ", restartMaxRestarts
  write (*, '(1X, A)') "    dt   iterations/solve (Krylov 30, bound; restarted, bound)   CN residual   " // &
    "unitarity   Cayley(-dt)Cayley(dt)-1   restarted - one cycle"
  do iDt = 1, nDt
    write (*, '(1X, F8.4, 4(1X, I6), 4(1X, ES12.3))') timeSteps(iDt), iterationsOneCycle(iDt), boundOneCycle(iDt), &
      iterationsRestart(iDt), boundRestart(iDt), residual(iDt), unitarityDev(iDt), roundTrip(iDt), agreement(iDt)
  end do

  call check(error, all(iterationsRestart > restartKrylovDim))
  if (allocated(error)) error stop "T_Method_FieldFreeCayleyFem failure: the restart path of the GMRES did not run"

  call check(error, all(residual < 1e-11_R64))
  if (allocated(error)) error stop "T_Method_FieldFreeCayleyFem failure: restarted Cayley step does not solve its system"

  call check(error, all(unitarityDev < 1e-11_R64))
  if (allocated(error)) error stop "T_Method_FieldFreeCayleyFem failure: restarted Cayley step not unitary"

  call check(error, all(roundTrip < 1e-11_R64))
  if (allocated(error)) error stop "T_Method_FieldFreeCayleyFem failure: Cayley(-dt) is not the inverse of Cayley(dt)"

  call check(error, all(agreement < 1e-11_R64))
  if (allocated(error)) error stop "T_Method_FieldFreeCayleyFem failure: restarted and one-cycle Cayley steps differ"

  call check(error, all(iterationsOneCycle <= boundOneCycle))
  if (allocated(error)) error stop "T_Method_FieldFreeCayleyFem failure: too many GMRES iterations (Krylov dimension 30)"

  call check(error, all(iterationsRestart <= boundRestart))
  if (allocated(error)) error stop "T_Method_FieldFreeCayleyFem failure: too many GMRES iterations (restarted)"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Two complex test orbitals at the nodes, normalized in the grid metric:
  !> exp(-r) exp(i (0.8 z - 0.3 y)) and (z + 0.5 x) exp(-r/2) exp(0.5 i x).
  subroutine FillTestOrbitals(orbs)
    complex(R64), intent(out), allocatable :: orbs(:, :)

    real(R64) :: x, y, z, r
    integer(I32) :: ix, iy, iz, node, k

    allocate (orbs(Grid_nPoints, N_ORBS))
    do iz = 1, Grid_Fem_nz
      do iy = 1, Grid_Fem_ny
        do ix = 1, Grid_Fem_nx
          node = Grid_Fem_Index(ix, iy, iz)
          x = Grid_Fem_axes(1) % points(ix)
          y = Grid_Fem_axes(2) % points(iy)
          z = Grid_Fem_axes(3) % points(iz)
          r = sqrt(x**2 + y**2 + z**2)
          orbs(node, 1) = exp(-r) * exp(cmplx(0.0_R64, 0.8_R64 * z - 0.3_R64 * y, kind=R64))
          orbs(node, 2) = (z + 0.5_R64 * x) * exp(-0.5_R64 * r) * exp(cmplx(0.0_R64, 0.5_R64 * x, kind=R64))
        end do
      end do
    end do
    do k = 1, N_ORBS
      orbs(:, k) = orbs(:, k) / MetricNorm(orbs(:, k))
    end do

  end subroutine

  !> One Cayley step of the orbitals with dt; nIterations = the largest GMRES
  !> iteration count of its solves.
  subroutine PropagateCounted(nIterations, orbs, dt)
    integer(I32), intent(out) :: nIterations
    complex(R64), intent(inout), contiguous :: orbs(:, :)
    real(R64), intent(in) :: dt

    Method_Mb_OrbBased_cayleyMaxIterations = 0
    call Method_Mb_OrbBased_PropagateFieldFree(orbs, dt, .false.)
    nIterations = Method_Mb_OrbBased_cayleyMaxIterations

  end subroutine

  !> |(1 + i dt h0/2) phiNew - (1 - i dt h0/2) phiOld| / |phiOld|.
  function CrankNicolsonResidual(phiNew, phiOld, dt) result(res)
    use M_Utils_Constants

    complex(R64), intent(in), contiguous :: phiNew(:)
    complex(R64), intent(in), contiguous :: phiOld(:)
    real(R64), intent(in) :: dt
    real(R64) :: res

    complex(R64), allocatable :: h0New(:), h0Old(:)

    allocate (h0New, h0Old, mold=phiNew)
    call ApplyH0(h0New, phiNew)
    call ApplyH0(h0Old, phiOld)
    res = MetricNorm(phiNew + 0.5_R64 * IU * dt * h0New - (phiOld - 0.5_R64 * IU * dt * h0Old)) / MetricNorm(phiOld)

  end function

  !> h0phi = (T + V_nuc) phi (body type 1; nodal potential of grid.fem).
  subroutine ApplyH0(h0phi, phiIn)
    complex(R64), intent(out), contiguous :: h0phi(:)
    complex(R64), intent(in), contiguous :: phiIn(:)

    call SysKinetic_MultiplyWithKineticOp(h0phi, phiIn, 0.0_R64, bt_=1)
    h0phi(:) = h0phi(:) + potential(:) * phiIn(:)

  end subroutine

  !> Overlap matrix <orbs_i|orbs_j> in the grid metric.
  function OverlapMatrix(orbs) result(res)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64) :: res(size(orbs, 2), size(orbs, 2))

    integer(I32) :: i, j

    do j = 1, size(orbs, 2)
      do i = 1, size(orbs, 2)
        res(i, j) = Grid_InnerProduct(orbs(:, i), orbs(:, j))
      end do
    end do

  end function

  !> Norm in the grid metric, sqrt(<f|f>).
  function MetricNorm(f) result(res)
    complex(R64), intent(in), contiguous :: f(:)
    real(R64) :: res

    res = sqrt(real(Grid_InnerProduct(f, f), kind=R64))

  end function

end program
