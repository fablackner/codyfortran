program T_Grid_Fem_Poisson
  !-----------------------------------------------------------------------------
  ! Galerkin Poisson solves on the FEM grid (mesh of T_Grid_Fem: box +-8,
  ! nLocals 6, centers (0, 0, -+0.7), hMin 0.2, growth 1.6, hMax 2;
  ! multipole boundary values up to l = 6 about the default center, the
  ! centroid of the centers = the origin; axis mix-ups, a non-default center
  ! and the l <= 12 boundary values are tested by T_Grid_Fem_Asymmetric):
  !  (a) unit Gaussian charge rho = pi^-1.5 exp(-r^2) at the origin: potential
  !      against erf(r)/r and the self-energy 1/2 sum load pot = 1/sqrt(2 pi)
  !  (b) the same charge centered at (0.5, -0.3, 0.8) (multipoles l > 0)
  !  (c) complex load: identical to real and imaginary parts solved separately
  !  (d) point charge q = -1 on the node (0, 0, 0.7): potential against -1/r,
  !      finite value at the nucleus
  !  (e) point charge q = 1 off the nodes: far field against 1/|r - R|
  ! Prototype references: (a) 5.903e-7 (r < 4), 0.398942323252; (b) 4.775e-6;
  ! (d) 6.669e-6 (r > 1.5), 1.269e-6 (r > 3), -229.5271624325 at the nucleus.
  !-----------------------------------------------------------------------------
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: nucleus(3) = [0.0_R64, 0.0_R64, 0.7_R64]
  real(R64), parameter :: offNode(3) = [0.13_R64, -0.21_R64, 0.37_R64]
  integer(I32), parameter :: nRepeat = 5

  integer(I32) :: ix, iy, iz, i, iRepeat
  integer(I64) :: tick0, tick1, tickRate
  real(R64) :: selfEnergyExact, errCentered, errOffCenter, errComplex, selfEnergy, potNucleus
  real(R64) :: errPoint15, errPoint3, errOffNode3, tPoisson, tPoissonComplex, tPointCharges
  real(R64), allocatable :: dist(:)
  complex(R64), allocatable :: load(:), pot(:), potRe(:), potIm(:)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Grid/Fem/T_Grid_Fem_Poisson.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Grid_Setup

  allocate (pot(Grid_nPoints), potRe(Grid_nPoints), potIm(Grid_nPoints))

  !==========================================
  call Say_Section("perform action")
  !==========================================

  !------------------------------------
  ! (a) unit Gaussian charge at the origin: V = erf(r)/r (2/sqrt(pi) at r = 0)
  !------------------------------------

  load = GaussianLoad([0.0_R64, 0.0_R64, 0.0_R64], 0.0_R64)
  call Grid_Fem_SolvePoisson(pot, load)

  dist = Distances([0.0_R64, 0.0_R64, 0.0_R64])
  errCentered = maxval(abs(real(pot, kind=R64) - ErfOverR(dist)), mask=dist < 4.0_R64)
  selfEnergy = 0.5_R64 * sum(real(load, kind=R64) * real(pot, kind=R64))

  call check(error, errCentered, 0.0_R64, thr=1.8e-6_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Poisson centered Gaussian potential failure"
  call check(error, selfEnergy, 1.0_R64 / sqrt(2.0_R64 * PI), thr=1.3e-7_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Poisson self-energy failure"
  call check(error, maxval(abs(aimag(pot))), 0.0_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Poisson real load with imaginary potential failure"

  !------------------------------------
  ! (b) off-center Gaussian charge (dipole and higher multipoles about the
  ! expansion center)
  !------------------------------------

  load = GaussianLoad([0.5_R64, -0.3_R64, 0.8_R64], 0.0_R64)
  call Grid_Fem_SolvePoisson(pot, load)

  dist = Distances([0.5_R64, -0.3_R64, 0.8_R64])
  errOffCenter = maxval(abs(real(pot, kind=R64) - ErfOverR(dist)))

  call check(error, errOffCenter, 0.0_R64, thr=1.0e-5_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Poisson off-center Gaussian potential failure"

  !------------------------------------
  ! (c) complex load (Gaussian times exp(0.3 i x)): the solve equals the real
  ! and imaginary parts solved separately
  !------------------------------------

  load = GaussianLoad([0.5_R64, -0.3_R64, 0.8_R64], 0.3_R64)
  call Grid_Fem_SolvePoisson(pot, load)
  call Grid_Fem_SolvePoisson(potRe, cmplx(real(load, kind=R64), 0.0_R64, kind=R64))
  call Grid_Fem_SolvePoisson(potIm, cmplx(aimag(load), 0.0_R64, kind=R64))

  errComplex = maxval(abs(pot - (potRe + (0.0_R64, 1.0_R64) * potIm))) / maxval(abs(pot))
  call check(error, errComplex, 0.0_R64, thr=3.0e-15_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Poisson complex load failure"

  !------------------------------------
  ! (d) point charge q = -1 on the node (0, 0, 0.7): V = -1/r away from it,
  ! finite at the nucleus node
  !------------------------------------

  call Grid_Fem_SolvePoissonPointCharges(pot, [-1.0_R64], reshape(nucleus, [3, 1]))

  dist = Distances(nucleus)
  errPoint15 = maxval(abs(real(pot, kind=R64) + 1.0_R64 / max(dist, 1.0e-12_R64)), mask=dist > 1.5_R64)
  errPoint3 = maxval(abs(real(pot, kind=R64) + 1.0_R64 / max(dist, 1.0e-12_R64)), mask=dist > 3.0_R64)

  potNucleus = 0.0_R64
  do iz = 1, Grid_Fem_nz
    do iy = 1, Grid_Fem_ny
      do ix = 1, Grid_Fem_nx
        i = Grid_Fem_Index(ix, iy, iz)
        if (dist(i) < tiny(1.0_R64)) potNucleus = real(pot(i), kind=R64)
      end do
    end do
  end do

  call check(error, errPoint15, 0.0_R64, thr=2.0e-5_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Poisson point charge (r > 1.5) failure"
  call check(error, errPoint3, 0.0_R64, thr=4.0e-6_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Poisson point charge (r > 3) failure"
  call check(error, potNucleus, -229.5271624325_R64, thr=1.0e-10_R64, rel=.true.)
  if (allocated(error)) error stop "T_Grid_Fem_Poisson point charge value at the nucleus failure"

  !------------------------------------
  ! (e) point charge q = 1 inside an element (Lagrange-interpolated load):
  ! far field 1/|r - R|
  !------------------------------------

  call Grid_Fem_SolvePoissonPointCharges(pot, [1.0_R64], reshape(offNode, [3, 1]))

  dist = Distances(offNode)
  errOffNode3 = maxval(abs(real(pot, kind=R64) - 1.0_R64 / max(dist, 1.0e-12_R64)), mask=dist > 3.0_R64)

  call check(error, errOffNode3, 0.0_R64, thr=3.0e-6_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Poisson off-node point charge failure"

  !------------------------------------
  ! timings (averages over nRepeat calls)
  !------------------------------------

  load = GaussianLoad([0.5_R64, -0.3_R64, 0.8_R64], 0.0_R64)
  call system_clock(tick0, tickRate)
  do iRepeat = 1, nRepeat
    call Grid_Fem_SolvePoisson(pot, load)
  end do
  call system_clock(tick1)
  tPoisson = real(tick1 - tick0, R64) / real(tickRate, R64) / nRepeat

  load = GaussianLoad([0.5_R64, -0.3_R64, 0.8_R64], 0.3_R64)
  call system_clock(tick0)
  do iRepeat = 1, nRepeat
    call Grid_Fem_SolvePoisson(pot, load)
  end do
  call system_clock(tick1)
  tPoissonComplex = real(tick1 - tick0, R64) / real(tickRate, R64) / nRepeat

  call system_clock(tick0)
  do iRepeat = 1, nRepeat
    call Grid_Fem_SolvePoissonPointCharges(pot, [-1.0_R64], reshape(nucleus, [3, 1]))
  end do
  call system_clock(tick1)
  tPointCharges = real(tick1 - tick0, R64) / real(tickRate, R64) / nRepeat

  !==========================================
  call Say_Section("printout results")
  !==========================================

  selfEnergyExact = 1.0_R64 / sqrt(2.0_R64 * PI)
  print *
  write (*, '(A, ES12.4)') " (a) centered Gaussian, max |V - erf(r)/r| (r < 4) = ", errCentered
  write (*, '(A, F16.12, A, ES11.3, A)') " (a) self-energy = ", selfEnergy, "  (error ", selfEnergy - selfEnergyExact, ")"
  write (*, '(A, ES12.4)') " (b) off-center Gaussian, max |V - erf(r)/r|       = ", errOffCenter
  write (*, '(A, ES12.4)') " (c) complex load vs separate parts, rel. diff.    = ", errComplex
  write (*, '(A, ES12.4)') " (d) point charge, max |V + 1/r| (r > 1.5)         = ", errPoint15
  write (*, '(A, ES12.4)') " (d) point charge, max |V + 1/r| (r > 3)           = ", errPoint3
  write (*, '(A, F16.10)') " (d) potential at the nucleus node                 = ", potNucleus
  write (*, '(A, ES12.4)') " (e) off-node charge, max |V - 1/r| (r > 3)        = ", errOffNode3
  write (*, '(A, F9.4, A)') " time SolvePoisson (real load)                     = ", tPoisson, " s"
  write (*, '(A, F9.4, A)') " time SolvePoisson (complex load)                  = ", tPoissonComplex, " s"
  write (*, '(A, F9.4, A)') " time SolvePoissonPointCharges                     = ", tPointCharges, " s"
  print *
  print *, "all fem Poisson checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Weighted load w rho of the unit Gaussian charge pi^-1.5 exp(-|r - r0|^2),
  !> times exp(i kappa x).
  function GaussianLoad(r0, kappa) result(field)
    real(R64), intent(in) :: r0(3)
    real(R64), intent(in) :: kappa
    complex(R64), allocatable :: field(:)

    integer(I32) :: jx, jy, jz, j
    real(R64) :: x, dr(3)

    allocate (field(Grid_nPoints))
    do jz = 1, Grid_Fem_nz
      do jy = 1, Grid_Fem_ny
        do jx = 1, Grid_Fem_nx
          j = Grid_Fem_Index(jx, jy, jz)
          x = Grid_Fem_axes(1) % points(jx)
          dr = [x, Grid_Fem_axes(2) % points(jy), Grid_Fem_axes(3) % points(jz)] - r0
          field(j) = Grid_Fem_weights(j) * PI**(-1.5_R64) * exp(-sum(dr**2)) * exp(cmplx(0.0_R64, kappa * x, kind=R64))
        end do
      end do
    end do

  end function

  !> Distances |r_i - r0| of all nodes.
  function Distances(r0) result(res)
    real(R64), intent(in) :: r0(3)
    real(R64), allocatable :: res(:)

    integer(I32) :: jx, jy, jz

    allocate (res(Grid_nPoints))
    do jz = 1, Grid_Fem_nz
      do jy = 1, Grid_Fem_ny
        do jx = 1, Grid_Fem_nx
          res(Grid_Fem_Index(jx, jy, jz)) = norm2([Grid_Fem_axes(1) % points(jx), Grid_Fem_axes(2) % points(jy), &
                                                   Grid_Fem_axes(3) % points(jz)] - r0)
        end do
      end do
    end do

  end function

  !> Potential erf(r)/r of the unit Gaussian charge (2/sqrt(pi) at r = 0).
  elemental function ErfOverR(r) result(v)
    real(R64), intent(in) :: r
    real(R64) :: v

    if (r > 1.0e-12_R64) then
      v = erf(r) / r
    else
      v = 2.0_R64 / sqrt(PI)
    end if

  end function

end program
