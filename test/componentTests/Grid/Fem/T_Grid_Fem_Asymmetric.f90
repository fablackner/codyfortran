program T_Grid_Fem_Asymmetric
  !-----------------------------------------------------------------------------
  ! FEM grid on a mesh whose three axes all differ, so that axis mix-ups show:
  ! nLocals 5, box [-7, 8.5] x [-6, 6.5] x [-8, 7], centers C1 = (0.4, -0.3,
  ! 0.9) (hMin 0.25) and C2 = (-0.8, 0.5, -0.6) (hMin 0.3), growth 1.6, hMax 2,
  ! graded x and z axes, explicit y bounds (first entry 1e-13 off the box,
  ! snapped), multipole boundary values up to l = 12 about (0.3, -0.2, 0.4)
  ! (neither the box center nor the centroid of the centers)
  ! -> n = (59, 47, 63), N = 174699.
  !  (a) mesh: node counts, graded x/z bounds of the prototype, explicit y
  !      bounds, centers at element vertices and nodes, convergence radius
  !  (b) quadrature and kinetic energy of a Gaussian
  !  (c) self-adjointness, eigenbasis round trip and diagonal form, shifted
  !      solves with a real and a complex shift
  !  (d) Gaussian charge at (0.5, -0.4, 0.2): potential and self-energy
  !  (e) boundary values: a node-delta load through Grid_Fem_SolvePoisson
  !      (multipoles l <= 12) against Grid_Fem_SolvePoissonPointCharges
  !      (exact boundary values) at the same node
  !  (f) point charge -1 on the C1 node: far field and value at the node
  !  (g) complex load: equal to real and imaginary parts solved separately
  !  (h) OpenMP: the operators called inside a parallel region equal serial
  !      calls; the k = 1 load is made wide on purpose (exp(-0.15 r^2), a
  !      fraction ~1e-2 of sum |load| beyond the convergence radius, above the
  !      warning threshold 1e-4): the divergence warning of
  !      Grid_Fem_SolvePoisson is expected, printed once from the parallel
  !      region (once-only flag in a critical section)
  !  (i) node coordinates (Grid_Fem_FillNodeCoordinate): the x, y and z unit
  !      vectors give the node coordinates exactly, a general direction d the
  !      projection d . r_i (layout x fastest)
  ! Prototype references (numpy, same mesh): (a) radius 5.800240710119;
  ! (b) 3.4626e-6, kinetic 1.500000952185; (d) 2.1195e-5, self-energy
  ! 0.398942065386; (e) 1.7988e-9; (f) 3.7289e-5 (r > 1.5), 4.2122e-6 (r > 3),
  ! -108.7928651665 at the node.
  !-----------------------------------------------------------------------------
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use testdrive, only: check, error_type
  !$ use omp_lib, only: omp_get_max_threads

  implicit none

  ! Element boundaries: graded x and z (prototype), explicit y (input, the
  ! first entry snapped onto the box)
  real(R64), parameter :: xBoundsRef(0:15) = [-7.0_R64, -5.061919485233116_R64, -3.160866422795566_R64, &
                                              -1.8995150898575044_R64, -1.1942917433383793_R64, -0.8_R64, &
                                              -0.4240780636800942_R64, 0.08673171829357808_R64, 0.4_R64, &
                                              0.6836213637513877_R64, 1.1603021601860586_R64, 1.9614556614766607_R64, &
                                              3.3079493982167723_R64, 5.0386252190781615_R64, 6.769312609539115_R64, 8.5_R64]
  real(R64), parameter :: yBoundsRef(0:12) = [-6.0_R64, -4.3_R64, -2.7_R64, -1.5_R64, -0.8_R64, -0.3_R64, 0.1_R64, &
                                              0.5_R64, 1.1_R64, 2.0_R64, 3.3_R64, 4.9_R64, 6.5_R64]
  real(R64), parameter :: zBoundsRef(0:16) = [-8.0_R64, -6.184932700889253_R64, -4.369865401778506_R64, &
                                              -2.661037139825657_R64, -1.585709303724614_R64, -0.9618903652812839_R64, &
                                              -0.6_R64, -0.2864697724447034_R64, 0.2136142835519304_R64, &
                                              0.6387248041474882_R64, 0.9_R64, 1.1742551031487782_R64, &
                                              1.6290280881217694_R64, 2.3831390240259136_R64, 3.633615332900102_R64, &
                                              5.31419851877206_R64, 7.0_R64]
  real(R64), parameter :: gaussCenter(3) = [0.5_R64, -0.4_R64, 0.2_R64]
  real(R64), parameter :: nucleus(3) = [0.4_R64, -0.3_R64, 0.9_R64]
  real(R64), parameter :: deltaTarget(3) = [1.2_R64, -1.0_R64, 1.6_R64]
  real(R64), parameter :: coordDirection(3) = [0.6_R64, -1.3_R64, 0.45_R64]
  real(R64), parameter :: shift = 0.7_R64
  complex(R64), parameter :: complexShift = (3.0_R64, -40.0_R64)
  integer(I32), parameter :: nFields = 4
  integer(I32), parameter :: nOps = 8

  integer(I32) :: ix, iy, iz, i, iAxis, iCenter, iDelta, k, iOp, nThreads, deltaIndices(3), nodeIndices(3)
  real(R64) :: center, integral, kinetic, relSymmetry, relRoundTrip, relDiagonal, relResidual, relResidualComplex
  real(R64) :: eigenvalue, diffMax, scaleMax, errGauss, selfEnergy, errDelta, errPoint15, errPoint3, potNucleus
  real(R64) :: errComplex, relParallel(nOps), deltaNode(3), unitVector(3), errCoordAxes, errCoordDirection
  complex(R64) :: lhs, rhs
  real(R64), allocatable :: dist(:), coord(:)
  complex(R64), allocatable :: g(:), lapG(:), f(:), h(:), lapF(:), lapH(:), c(:), back(:), cLap(:), u(:), lapU(:)
  complex(R64), allocatable :: load(:), pot(:), potRe(:), potIm(:), fields(:, :), outParallel(:, :, :), work(:)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Grid/Fem/T_Grid_Fem_Asymmetric.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Grid_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  !------------------------------------
  ! (a) mesh: three different axes, element boundaries of the prototype (x,
  ! z) and of the input (y), every center coordinate an element boundary and
  ! a node, the expansion center as given
  !------------------------------------

  call check(error, Grid_Fem_nx, 59)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric nx failure"
  call check(error, Grid_Fem_ny, 47)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric ny failure"
  call check(error, Grid_Fem_nz, 63)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric nz failure"
  call check(error, Grid_nPoints, 174699)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric Grid_nPoints failure"
  call check(error, Grid_Fem_Index(Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz), Grid_nPoints)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric Grid_Fem_Index failure"

  call check(error, size(Grid_Fem_axes(1) % bounds), size(xBoundsRef))
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric x element count failure"
  call check(error, maxval(abs(Grid_Fem_axes(1) % bounds - xBoundsRef)), 0.0_R64, thr=1.0e-12_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric x element boundaries failure"
  call check(error, size(Grid_Fem_axes(2) % bounds), size(yBoundsRef))
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric y element count failure"
  call check(error, maxval(abs(Grid_Fem_axes(2) % bounds - yBoundsRef)) < tiny(1.0_R64), &
             "explicit y bounds differ from the input")
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric explicit y bounds failure"
  call check(error, size(Grid_Fem_axes(3) % bounds), size(zBoundsRef))
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric z element count failure"
  call check(error, maxval(abs(Grid_Fem_axes(3) % bounds - zBoundsRef)), 0.0_R64, thr=1.0e-12_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric z element boundaries failure"

  do iCenter = 1, Grid_Fem_nCenters
    do iAxis = 1, 3
      center = Grid_Fem_centers(iAxis, iCenter)
      call check(error, minval(abs(Grid_Fem_axes(iAxis) % bounds - center)) < tiny(1.0_R64), &
                 "center coordinate is not an element boundary")
      if (allocated(error)) error stop "T_Grid_Fem_Asymmetric center bound failure"
      call check(error, minval(abs(Grid_Fem_axes(iAxis) % points - center)) < tiny(1.0_R64), &
                 "center coordinate is not a node")
      if (allocated(error)) error stop "T_Grid_Fem_Asymmetric center node failure"
    end do
  end do

  call check(error, maxval(abs(Grid_Fem_multipoleCenter - [0.3_R64, -0.2_R64, 0.4_R64])) < tiny(1.0_R64), &
             "multipole center differs from the input")
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric multipole center failure"
  call check(error, Grid_Fem_multipoleRadius, 5.800240710119_R64, thr=1.0e-11_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric multipole convergence radius failure"

  !------------------------------------
  ! (b) quadrature and kinetic energy of g = exp(-|r - r0|^2):
  ! sum w g = pi^1.5, -1/2 <g|lap g>/<g|g> = 3/2
  !------------------------------------

  g = Gaussian(1.0_R64, gaussCenter, 0.0_R64)

  integral = sum(Grid_Fem_weights * real(g, kind=R64))
  call check(error, integral, PI**1.5_R64, thr=1.0e-5_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric Gaussian quadrature failure"

  allocate (lapG(Grid_nPoints))
  call Grid_Fem_ApplyLaplacian(lapG, g)

  kinetic = -0.5_R64 * real(Grid_InnerProduct(g, lapG), kind=R64) / real(Grid_InnerProduct(g, g), kind=R64)
  call check(error, kinetic, 1.5_R64, thr=3.0e-6_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric kinetic energy failure"
  call check(error, kinetic, 1.500000952185_R64, thr=1.0e-10_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric kinetic energy vs prototype failure"

  !------------------------------------
  ! (c) self-adjointness, eigenbasis and shifted solves for complex
  ! off-center functions
  !------------------------------------

  f = Gaussian(0.8_R64, [0.3_R64, -0.2_R64, 0.5_R64], 0.3_R64)
  h = Gaussian(1.1_R64, [-0.4_R64, 0.1_R64, -0.6_R64], -0.5_R64)

  allocate (lapF(Grid_nPoints), lapH(Grid_nPoints))
  call Grid_Fem_ApplyLaplacian(lapF, f)
  call Grid_Fem_ApplyLaplacian(lapH, h)

  lhs = Grid_InnerProduct(f, lapH)
  rhs = Grid_InnerProduct(lapF, h)
  relSymmetry = abs(lhs - rhs) / abs(lhs)
  call check(error, relSymmetry, 0.0_R64, thr=1.0e-13_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric Laplacian symmetry failure"

  allocate (c(Grid_nPoints), back(Grid_nPoints), cLap(Grid_nPoints))
  call Grid_Fem_ToEigenbasis(c, f)
  call Grid_Fem_FromEigenbasis(back, c)

  relRoundTrip = maxval(abs(back - f)) / maxval(abs(f))
  call check(error, relRoundTrip, 0.0_R64, thr=2.0e-14_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric eigenbasis round trip failure"

  ! ToEigenbasis(lap f) = -D ToEigenbasis(f), D(i, j, k) = lx(i) + ly(j) + lz(k)
  call Grid_Fem_ToEigenbasis(cLap, lapF)

  diffMax = 0.0_R64
  scaleMax = 0.0_R64
  do iz = 1, Grid_Fem_nz
    do iy = 1, Grid_Fem_ny
      do ix = 1, Grid_Fem_nx
        i = Grid_Fem_Index(ix, iy, iz)
        eigenvalue = Grid_Fem_axes(1) % eigenvalues(ix) + Grid_Fem_axes(2) % eigenvalues(iy) + &
                     Grid_Fem_axes(3) % eigenvalues(iz)
        diffMax = max(diffMax, abs(cLap(i) + eigenvalue * c(i)))
        scaleMax = max(scaleMax, abs(eigenvalue * c(i)))
      end do
    end do
  end do

  relDiagonal = diffMax / scaleMax
  call check(error, relDiagonal, 0.0_R64, thr=2.0e-12_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric eigenbasis diagonalization failure"

  allocate (u(Grid_nPoints), lapU(Grid_nPoints))
  call Grid_Fem_SolveShiftedLaplacian(u, f, shift)
  call Grid_Fem_ApplyLaplacian(lapU, u)

  relResidual = maxval(abs(-lapU + shift * u - f)) / maxval(abs(f))
  call check(error, relResidual, 0.0_R64, thr=5.0e-12_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric shifted solve failure"

  call Grid_Fem_SolveShiftedLaplacian(u, f, complexShift)
  call Grid_Fem_ApplyLaplacian(lapU, u)

  relResidualComplex = maxval(abs(-lapU + complexShift * u - f)) / maxval(abs(f))
  call check(error, relResidualComplex, 0.0_R64, thr=5.0e-13_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric complex-shift solve failure"

  !------------------------------------
  ! (d) unit Gaussian charge pi^-1.5 exp(-|r - r0|^2): V = erf(r)/r,
  ! self-energy 1/sqrt(2 pi)
  !------------------------------------

  allocate (pot(Grid_nPoints), potRe(Grid_nPoints), potIm(Grid_nPoints))

  load = Grid_Fem_weights * PI**(-1.5_R64) * g
  call Grid_Fem_SolvePoisson(pot, load)

  dist = Distances(gaussCenter)
  errGauss = maxval(abs(real(pot, kind=R64) - ErfOverR(dist)))
  selfEnergy = 0.5_R64 * sum(real(load, kind=R64) * real(pot, kind=R64))

  call check(error, errGauss, 0.0_R64, thr=6.0e-5_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric Gaussian potential failure"
  call check(error, selfEnergy, 1.0_R64 / sqrt(2.0_R64 * PI), thr=6.5e-7_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric self-energy failure"
  call check(error, selfEnergy, 0.398942065386_R64, thr=1.0e-11_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric self-energy vs prototype failure"

  !------------------------------------
  ! (e) boundary values: node-delta load at the node nearest to (1.2, -1.0,
  ! 1.6) (|r - c| = 1.66 of the convergence radius 5.80) through the
  ! multipole path against the exact boundary values of the point-charge
  ! path (same Galerkin load 4 pi delta): only V_B differs
  !------------------------------------

  do iAxis = 1, 3
    deltaIndices(iAxis) = minloc(abs(Grid_Fem_axes(iAxis) % points - deltaTarget(iAxis)), 1)
    deltaNode(iAxis) = Grid_Fem_axes(iAxis) % points(deltaIndices(iAxis))
  end do
  iDelta = Grid_Fem_Index(deltaIndices(1), deltaIndices(2), deltaIndices(3))

  load = (0.0_R64, 0.0_R64)
  load(iDelta) = (1.0_R64, 0.0_R64)
  call Grid_Fem_SolvePoisson(potRe, load)
  call Grid_Fem_SolvePoissonPointCharges(potIm, [1.0_R64], reshape(deltaNode, [3, 1]))

  errDelta = maxval(abs(potRe - potIm))
  call check(error, errDelta, 0.0_R64, thr=5.0e-9_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric multipole boundary values failure"

  !------------------------------------
  ! (f) point charge q = -1 on the C1 node: V = -1/r away from it, finite at
  ! the node
  !------------------------------------

  call Grid_Fem_SolvePoissonPointCharges(pot, [-1.0_R64], reshape(nucleus, [3, 1]))

  dist = Distances(nucleus)
  errPoint15 = maxval(abs(real(pot, kind=R64) + 1.0_R64 / max(dist, 1.0e-12_R64)), mask=dist > 1.5_R64)
  errPoint3 = maxval(abs(real(pot, kind=R64) + 1.0_R64 / max(dist, 1.0e-12_R64)), mask=dist > 3.0_R64)
  potNucleus = real(pot(minloc(dist, 1)), kind=R64)

  call check(error, minval(dist) < tiny(1.0_R64), "C1 is not a node")
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric C1 node failure"
  call check(error, errPoint15, 0.0_R64, thr=1.1e-4_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric point charge (r > 1.5) failure"
  call check(error, errPoint3, 0.0_R64, thr=1.3e-5_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric point charge (r > 3) failure"
  call check(error, potNucleus, -108.7928651665_R64, thr=1.0e-10_R64, rel=.true.)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric point charge value at the nucleus failure"

  !------------------------------------
  ! (g) complex load (Gaussian charge times exp(0.3 i x)): equal to the real
  ! and imaginary parts solved separately
  !------------------------------------

  load = Grid_Fem_weights * PI**(-1.5_R64) * Gaussian(1.0_R64, gaussCenter, 0.3_R64)
  call Grid_Fem_SolvePoisson(pot, load)
  call Grid_Fem_SolvePoisson(potRe, cmplx(real(load, kind=R64), 0.0_R64, kind=R64))
  call Grid_Fem_SolvePoisson(potIm, cmplx(aimag(load), 0.0_R64, kind=R64))

  errComplex = maxval(abs(pot - (potRe + (0.0_R64, 1.0_R64) * potIm))) / maxval(abs(pot))
  call check(error, errComplex, 0.0_R64, thr=3.0e-15_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric complex load failure"

  !------------------------------------
  ! (h) OpenMP: every operator called for nFields inputs inside a parallel
  ! region (where their own OpenMP loops are off) equals the serial call
  !------------------------------------

  nThreads = 1
  !$ nThreads = omp_get_max_threads()

  allocate (fields(Grid_nPoints, nFields), outParallel(Grid_nPoints, nOps, nFields), work(Grid_nPoints))
  do k = 1, nFields
    fields(:, k) = Gaussian(0.5_R64 + 0.1_R64 * k, [0.1_R64 * k, -0.2_R64, 0.15_R64 * k], 0.3_R64 * (k - 2))
  end do
  ! wide on purpose: raises the once-only multipole warning inside the region
  fields(:, 1) = Gaussian(0.15_R64, [0.1_R64, -0.2_R64, 0.15_R64], -0.3_R64)

  !$omp parallel do default(shared) private(k, iOp) schedule(static, 1)
  do k = 1, nFields
    do iOp = 1, nOps
      call ApplyOperator(outParallel(:, iOp, k), fields(:, k), iOp, k)
    end do
  end do
  !$omp end parallel do

  relParallel = 0.0_R64
  do k = 1, nFields
    do iOp = 1, nOps
      call ApplyOperator(work, fields(:, k), iOp, k)
      relParallel(iOp) = max(relParallel(iOp), maxval(abs(outParallel(:, iOp, k) - work)) / maxval(abs(work)))
    end do
  end do

  call check(error, maxval(relParallel), 0.0_R64, thr=5.0e-15_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric parallel vs serial calls failure"

  !------------------------------------
  ! (i) node coordinates: exact along the three axes (three different
  ! meshes, so an axis mix-up shows), d . r_i along a general direction
  !------------------------------------

  errCoordAxes = 0.0_R64
  do iAxis = 1, 3
    unitVector(:) = 0.0_R64
    unitVector(iAxis) = 1.0_R64
    call Grid_Fem_FillNodeCoordinate(coord, unitVector)
    call check(error, size(coord), Grid_nPoints)
    if (allocated(error)) error stop "T_Grid_Fem_Asymmetric node coordinate size failure"
    do iz = 1, Grid_Fem_nz
      do iy = 1, Grid_Fem_ny
        do ix = 1, Grid_Fem_nx
          nodeIndices = [ix, iy, iz]
          errCoordAxes = max(errCoordAxes, abs(coord(Grid_Fem_Index(ix, iy, iz)) - &
                                               Grid_Fem_axes(iAxis) % points(nodeIndices(iAxis))))
        end do
      end do
    end do
  end do

  call Grid_Fem_FillNodeCoordinate(coord, coordDirection)
  errCoordDirection = 0.0_R64
  do iz = 1, Grid_Fem_nz
    do iy = 1, Grid_Fem_ny
      do ix = 1, Grid_Fem_nx
        errCoordDirection = max(errCoordDirection, abs(coord(Grid_Fem_Index(ix, iy, iz)) - &
                                dot_product(coordDirection, [Grid_Fem_axes(1) % points(ix), &
                                                             Grid_Fem_axes(2) % points(iy), &
                                                             Grid_Fem_axes(3) % points(iz)])))
      end do
    end do
  end do

  call check(error, errCoordAxes < tiny(1.0_R64), "node coordinates along the axes must be exact")
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric node coordinates (axes) failure"
  call check(error, errCoordDirection, 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Asymmetric node coordinates (direction) failure"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(A, 3I4, A, I0)') " (a) n = ", Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, ", N = ", Grid_nPoints
  write (*, '(A, F16.12)') " (a) multipole convergence radius              = ", Grid_Fem_multipoleRadius
  write (*, '(A, ES12.4)') " (b) sum w g - pi^1.5                          = ", integral - PI**1.5_R64
  write (*, '(A, F16.12)') " (b) kinetic energy of g                       = ", kinetic
  write (*, '(A, ES12.4)') " (c) |<f|lap h> - <lap f|h>| / |<f|lap h>|      = ", relSymmetry
  write (*, '(A, ES12.4)') " (c) round trip, max rel. error                = ", relRoundTrip
  write (*, '(A, ES12.4)') " (c) ToEigenbasis(lap f) + D c, rel. error     = ", relDiagonal
  write (*, '(A, ES12.4)') " (c) shifted solve (0.7), rel. residual        = ", relResidual
  write (*, '(A, ES12.4)') " (c) shifted solve (3 - 40i), rel. residual    = ", relResidualComplex
  write (*, '(A, ES12.4)') " (d) Gaussian charge, max |V - erf(r)/r|       = ", errGauss
  write (*, '(A, F16.12, A, ES11.3, A)') " (d) self-energy = ", selfEnergy, "  (error ", &
    selfEnergy - 1.0_R64 / sqrt(2.0_R64 * PI), ")"
  write (*, '(A, 3F9.4)') " (e) node-delta load at ", deltaNode
  write (*, '(A, ES12.4)') " (e) max |V(multipole V_B) - V(exact V_B)|     = ", errDelta
  write (*, '(A, ES12.4)') " (f) point charge, max |V + 1/r| (r > 1.5)     = ", errPoint15
  write (*, '(A, ES12.4)') " (f) point charge, max |V + 1/r| (r > 3)       = ", errPoint3
  write (*, '(A, F16.10)') " (f) potential at the nucleus node             = ", potNucleus
  write (*, '(A, ES12.4)') " (g) complex load vs separate parts, rel. diff = ", errComplex
  write (*, '(A, I0, A)') " (h) parallel vs serial calls (", nThreads, " threads), max rel. diff per operator:"
  write (*, '(A, 8ES9.1)') "     ", relParallel
  if (nThreads < 2) print *, "    (one OpenMP thread only: the parallel region ran serially)"
  print *, "    (the multipole divergence warning above is expected: the wide k = 1 load of (h) reaches beyond the radius)"
  write (*, '(A, 2ES12.4)') " (i) node coordinates, max error axes / d . r  = ", errCoordAxes, errCoordDirection
  print *
  print *, "all fem asymmetric-mesh checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Operator iOp applied to the input f of field k (see (h)).
  subroutine ApplyOperator(res, f, iOp, k)
    complex(R64), intent(out), contiguous :: res(:)
    complex(R64), intent(in), contiguous  :: f(:)
    integer(I32), intent(in)              :: iOp
    integer(I32), intent(in)              :: k

    select case (iOp)
    case (1)
      call Grid_Fem_ApplyLaplacian(res, f)
    case (2)
      call Grid_Fem_ToEigenbasis(res, f)
    case (3)
      call Grid_Fem_FromEigenbasis(res, f)
    case (4)
      call Grid_Fem_SolveShiftedLaplacian(res, f, 0.3_R64 * k)
    case (5)
      call Grid_Fem_SolveShiftedLaplacian(res, f, cmplx(1.0_R64 * k, -10.0_R64 * k, kind=R64))
    case (6)
      call Grid_Fem_SolvePoisson(res, f * Grid_Fem_weights)
    case (7)
      call Grid_Fem_SolvePoisson(res, cmplx(real(f, kind=R64) * Grid_Fem_weights, 0.0_R64, kind=R64))
    case default
      call Grid_Fem_SolvePoissonPointCharges(res, [-1.0_R64 * k, 0.5_R64], &
                                             reshape([0.11_R64 * k, -0.2_R64, 0.3_R64, 1.0_R64, 0.0_R64, -0.37_R64 * k], &
                                                     [3, 2]))
    end select

  end subroutine

  !> exp(-alpha |r - r0|^2) exp(i kappa (x - x0)) on the nodes.
  function Gaussian(alpha, r0, kappa) result(field)
    real(R64), intent(in) :: alpha
    real(R64), intent(in) :: r0(3)
    real(R64), intent(in) :: kappa
    complex(R64), allocatable :: field(:)

    integer(I32) :: jx, jy, jz
    real(R64) :: d(3)

    allocate (field(Grid_nPoints))
    do jz = 1, Grid_Fem_nz
      do jy = 1, Grid_Fem_ny
        do jx = 1, Grid_Fem_nx
          d = [Grid_Fem_axes(1) % points(jx), Grid_Fem_axes(2) % points(jy), Grid_Fem_axes(3) % points(jz)] - r0
          field(Grid_Fem_Index(jx, jy, jz)) = exp(-alpha * sum(d**2)) * exp(cmplx(0.0_R64, kappa * d(1), kind=R64))
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
