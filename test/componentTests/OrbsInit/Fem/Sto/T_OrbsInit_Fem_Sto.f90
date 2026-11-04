!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Component test for the Slater-type-orbital initializer orbsInit.fem.sto on
!> a two-center mesh (A = (0, 0, 0), B = (0, 0, 1.4); box +-8, nLocals 6,
!> hMin 0.3, growth 1.6, hMax 2). Orbitals (unnormalized terms):
!>   1: exp(-rA)                                  1s on A
!>   2: exp(-1.2 rA) - 0.6 exp(-0.9 rB)           LCAO of two 1s
!>   3: x exp(-0.8 rA)                            2p_x (S_11 = x)
!>   4: sqrt(3) x y exp(-1.1 rA)                  3d_xy (S_2-2)
!>   5: rB^2 exp(-0.7 rB)                         3s on B (rho^(n-l-1))
!>   6: sqrt(15)/2 (x^2 - y^2) z exp(-1.3 rA)     4f (S_32)
!>
!> Checks:
!> - six orbitals defined, each normalized in the grid metric
!> - every orbital equals the analytic function normalized on the grid
!> - 1s at the nucleus node against the analytically normalized 1/sqrt(pi)
!> - LCAO sum (orbital 2) and the p, d and f reflection symmetries about A
!> - the same orbitals for every body type (and without bt_)
program T_OrbsInit_Fem_Sto
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use M_OrbsInit
  use M_OrbsInit_Fem
  use testdrive, only: check, error_type

  implicit none

  integer(I32), parameter :: nOrbs = 6
  real(R64), parameter :: centerB(3) = [0.0_R64, 0.0_R64, 1.4_R64]

  integer(I32) :: ix, iy, iz, i, iOrb, nx, ny, nodeA
  real(R64) :: x, y, z, rA, rB, errNorm, errAnalytic(nOrbs), errNucleus, errBodyType, errAxes
  real(R64) :: errPx, errPy, errDxy, errDyx, errDswap, errFswap, errFx, errSigma
  complex(R64), allocatable :: orbs(:, :), expected(:, :), orbBt2(:), orbNoBt(:)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/OrbsInit/Fem/Sto/T_OrbsInit_Fem_Sto.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call OrbsInit_Fabricate

  call Grid_Setup
  call OrbsInit_Setup

  nx = Grid_Fem_nx
  ny = Grid_Fem_ny

  allocate (orbs(Grid_nPoints, nOrbs), expected(Grid_nPoints, nOrbs), orbBt2(Grid_nPoints), orbNoBt(Grid_nPoints))

  !==========================================
  call Say_Section("perform action")
  !==========================================

  call check(error, OrbsInit_Fem_nOrbitals, nOrbs)
  if (allocated(error)) error stop "T_OrbsInit_Fem_Sto number of orbitals failure"

  do iOrb = 1, nOrbs
    call OrbsInit_InitializeOrb(orbs(:, iOrb), iOrb, bt_=1)
  end do

  !------------------------------------
  ! normalization
  !------------------------------------

  errNorm = 0.0_R64
  do iOrb = 1, nOrbs
    errNorm = max(errNorm, abs(Grid_InnerProduct(orbs(:, iOrb), orbs(:, iOrb)) - 1.0_R64))
  end do

  call check(error, errNorm, 0.0_R64, thr=1.0e-12_R64)
  if (allocated(error)) error stop "T_OrbsInit_Fem_Sto normalization failure"

  !------------------------------------
  ! analytic functions, normalized on the grid
  !------------------------------------

  do iz = 1, Grid_Fem_nz
    do iy = 1, ny
      do ix = 1, nx
        i = Grid_Fem_Index(ix, iy, iz)
        x = Grid_Fem_axes(1) % points(ix)
        y = Grid_Fem_axes(2) % points(iy)
        z = Grid_Fem_axes(3) % points(iz)
        rA = sqrt(x * x + y * y + z * z)
        rB = norm2([x, y, z] - centerB)

        expected(i, 1) = exp(-rA)
        expected(i, 2) = exp(-1.2_R64 * rA) - 0.6_R64 * exp(-0.9_R64 * rB)
        expected(i, 3) = x * exp(-0.8_R64 * rA)
        expected(i, 4) = sqrt(3.0_R64) * x * y * exp(-1.1_R64 * rA)
        expected(i, 5) = rB**2 * exp(-0.7_R64 * rB)
        expected(i, 6) = 0.5_R64 * sqrt(15.0_R64) * (x * x - y * y) * z * exp(-1.3_R64 * rA)
      end do
    end do
  end do

  do iOrb = 1, nOrbs
    expected(:, iOrb) = expected(:, iOrb) / sqrt(real(Grid_InnerProduct(expected(:, iOrb), expected(:, iOrb)), kind=R64))
    errAnalytic(iOrb) = maxval(abs(orbs(:, iOrb) - expected(:, iOrb))) / maxval(abs(expected(:, iOrb)))
  end do

  call check(error, maxval(errAnalytic), 0.0_R64, thr=1.0e-13_R64)
  if (allocated(error)) error stop "T_OrbsInit_Fem_Sto analytic orbitals failure"

  ! 1s at the nucleus node against the analytic normalization (zeta^3/pi)^(1/2)
  nodeA = Grid_Fem_Index(minloc(abs(Grid_Fem_axes(1) % points), dim=1), minloc(abs(Grid_Fem_axes(2) % points), dim=1), &
                         minloc(abs(Grid_Fem_axes(3) % points), dim=1))
  errNucleus = abs(real(orbs(nodeA, 1), kind=R64) - 1.0_R64 / sqrt(PI)) * sqrt(PI)

  call check(error, errNucleus, 0.0_R64, thr=1.0e-5_R64)
  if (allocated(error)) error stop "T_OrbsInit_Fem_Sto 1s value at the nucleus failure"

  !------------------------------------
  ! reflection symmetries about A (the x and y axes are mirror-symmetric)
  !------------------------------------

  errAxes = max(maxval(abs(Grid_Fem_axes(1) % points(nx:1:-1) + Grid_Fem_axes(1) % points)), &
                maxval(abs(Grid_Fem_axes(2) % points - Grid_Fem_axes(1) % points)))

  call check(error, nx .eq. ny .and. errAxes < 1.0e-12_R64, "the x and y axes must be identical and mirror-symmetric")
  if (allocated(error)) error stop "T_OrbsInit_Fem_Sto mesh symmetry failure"

  errSigma = 0.0_R64
  errPx = 0.0_R64
  errPy = 0.0_R64
  errDxy = 0.0_R64
  errDyx = 0.0_R64
  errDswap = 0.0_R64
  errFswap = 0.0_R64
  errFx = 0.0_R64
  do iz = 1, Grid_Fem_nz
    do iy = 1, ny
      do ix = 1, nx
        i = Grid_Fem_Index(ix, iy, iz)
        ! sigma-type LCAO and p_x under x -> -x and y -> -y
        errSigma = max(errSigma, abs(orbs(Grid_Fem_Index(nx + 1 - ix, iy, iz), 2) - orbs(i, 2)))
        errPx = max(errPx, abs(orbs(Grid_Fem_Index(nx + 1 - ix, iy, iz), 3) + orbs(i, 3)))
        errPy = max(errPy, abs(orbs(Grid_Fem_Index(ix, ny + 1 - iy, iz), 3) - orbs(i, 3)))
        ! d_xy: odd in x and in y, even under x <-> y
        errDxy = max(errDxy, abs(orbs(Grid_Fem_Index(nx + 1 - ix, iy, iz), 4) + orbs(i, 4)))
        errDyx = max(errDyx, abs(orbs(Grid_Fem_Index(ix, ny + 1 - iy, iz), 4) + orbs(i, 4)))
        errDswap = max(errDswap, abs(orbs(Grid_Fem_Index(iy, ix, iz), 4) - orbs(i, 4)))
        ! f (x^2 - y^2) z: odd under x <-> y, even in x
        errFswap = max(errFswap, abs(orbs(Grid_Fem_Index(iy, ix, iz), 6) + orbs(i, 6)))
        errFx = max(errFx, abs(orbs(Grid_Fem_Index(nx + 1 - ix, iy, iz), 6) - orbs(i, 6)))
      end do
    end do
  end do

  errSigma = errSigma / maxval(abs(orbs(:, 2)))
  errPx = errPx / maxval(abs(orbs(:, 3)))
  errPy = errPy / maxval(abs(orbs(:, 3)))
  errDxy = errDxy / maxval(abs(orbs(:, 4)))
  errDyx = errDyx / maxval(abs(orbs(:, 4)))
  errDswap = errDswap / maxval(abs(orbs(:, 4)))
  errFswap = errFswap / maxval(abs(orbs(:, 6)))
  errFx = errFx / maxval(abs(orbs(:, 6)))

  call check(error, max(errSigma, errPx, errPy), 0.0_R64, thr=1.0e-12_R64)
  if (allocated(error)) error stop "T_OrbsInit_Fem_Sto sigma and p symmetry failure"
  call check(error, max(errDxy, errDyx, errDswap), 0.0_R64, thr=1.0e-12_R64)
  if (allocated(error)) error stop "T_OrbsInit_Fem_Sto d symmetry failure"
  call check(error, max(errFswap, errFx), 0.0_R64, thr=1.0e-12_R64)
  if (allocated(error)) error stop "T_OrbsInit_Fem_Sto f symmetry failure"

  !------------------------------------
  ! the same orbitals for every body type and without bt_ (method.sb)
  !------------------------------------

  call OrbsInit_InitializeOrb(orbBt2, 4, bt_=2)
  call OrbsInit_InitializeOrb(orbNoBt, 4)
  errBodyType = max(maxval(abs(orbBt2 - orbs(:, 4))), maxval(abs(orbNoBt - orbs(:, 4))))

  call check(error, errBodyType <= 0.0_R64, "the orbitals must not depend on the body type")
  if (allocated(error)) error stop "T_OrbsInit_Fem_Sto body type failure"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(A, I0)') " Grid_nPoints                                   = ", Grid_nPoints
  write (*, '(A, ES12.4)') " max |<phi|phi> - 1|                            = ", errNorm
  write (*, '(A, 6ES10.2)') " max|phi - analytic| / max|phi|, orbitals 1-6   = ", errAnalytic
  write (*, '(A, ES12.4)') " 1s at the nucleus vs 1/sqrt(pi), rel.          = ", errNucleus
  write (*, '(A, 3ES12.4)') " sigma (x), p_x (x, y) reflections, rel.        = ", errSigma, errPx, errPy
  write (*, '(A, 3ES12.4)') " d_xy (x, y, x<->y) reflections, rel.           = ", errDxy, errDyx, errDswap
  write (*, '(A, 2ES12.4)') " f (x<->y, x) reflections, rel.                 = ", errFswap, errFx
  print *
  print *, "all fem STO initialization checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
