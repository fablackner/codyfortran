program T_Grid_Fem
  !-----------------------------------------------------------------------------
  ! FEM grid (grid.fem) on the graded H2 mesh of the numpy prototype: box +-8,
  ! nLocals 6, centers (0, 0, -0.7) and (0, 0, 0.7), hMin 0.2, growth 1.6,
  ! hMax 2 (n = (69, 69, 79), N = 376119).
  !  (a) mesh: node counts, element boundaries, centers on bounds and nodes
  !  (b) quadrature: sum w exp(-r^2) = pi^1.5
  !  (c) kinetic energy of exp(-r^2): -1/2 <g|lap g>/<g|g> = 3/2
  !  (d) self-adjointness: <f|lap g> = <lap f|g>
  !  (e) eigenbasis: round trip and -lap = S3 D S3^T W3
  !  (f) shifted solve: (-lap + s) u = f, real s and complex s (generic)
  !  (g) Grid_Orthonormalize with the FEM metric
  ! Reference numbers (prototype): (b) error 1.1e-8, (c) 1.4999993337.
  !-----------------------------------------------------------------------------
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use testdrive, only: check, error_type

  implicit none

  ! Element boundaries of the prototype mesh, rounded to 4 decimals
  real(R64), parameter :: xBoundsRef(0:14) = [-8.0_R64, -6.1892_R64, -4.3785_R64, -2.5946_R64, -1.3674_R64, &
                                              -0.6546_R64, -0.2405_R64, 0.0_R64, 0.2405_R64, 0.6546_R64, &
                                              1.3674_R64, 2.5946_R64, 4.3785_R64, 6.1892_R64, 8.0_R64]
  real(R64), parameter :: zBoundsRef(0:16) = [-8.0_R64, -6.0041_R64, -4.0082_R64, -2.3758_R64, -1.4706_R64, &
                                              -0.9733_R64, -0.7_R64, -0.4464_R64, 0.0_R64, 0.4464_R64, 0.7_R64, &
                                              0.9733_R64, 1.4706_R64, 2.3758_R64, 4.0082_R64, 6.0041_R64, 8.0_R64]
  real(R64), parameter :: shift = 0.7_R64
  complex(R64), parameter :: complexShift = (3.0_R64, -40.0_R64)
  integer(I32), parameter :: nRepeat = 5

  integer(I32) :: ix, iy, iz, i, iAxis, iCenter, iRepeat
  integer(I64) :: tick0, tick1, tickRate
  real(R64) :: center, integral, kinetic, relSymmetry, relRoundTrip, relDiagonal, relResidual, relResidualComplex
  real(R64) :: eigenvalue, diffMax, scaleMax, tLaplacian, tToEigen, tShifted
  complex(R64) :: lhs, rhs, ovlp12, norm11, norm22
  complex(R64), allocatable :: g(:), lapG(:), f(:), h(:), lapF(:), lapH(:), c(:), back(:), cLap(:)
  complex(R64), allocatable :: u(:), lapU(:), orbSet(:, :)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Grid/Fem/T_Grid_Fem.json", relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Grid_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  !------------------------------------
  ! (a) mesh: node counts and layout, element boundaries of the prototype,
  ! every center coordinate an element boundary and a node (exactly)
  !------------------------------------

  call check(error, Grid_Fem_nx, 69)
  if (allocated(error)) error stop "T_Grid_Fem nx failure"
  call check(error, Grid_Fem_ny, 69)
  if (allocated(error)) error stop "T_Grid_Fem ny failure"
  call check(error, Grid_Fem_nz, 79)
  if (allocated(error)) error stop "T_Grid_Fem nz failure"
  call check(error, Grid_nPoints, 376119)
  if (allocated(error)) error stop "T_Grid_Fem Grid_nPoints failure"
  call check(error, Grid_Fem_Index(Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz), Grid_nPoints)
  if (allocated(error)) error stop "T_Grid_Fem Grid_Fem_Index failure"

  call check(error, size(Grid_Fem_axes(1) % bounds), size(xBoundsRef))
  if (allocated(error)) error stop "T_Grid_Fem x element count failure"
  call check(error, maxval(abs(Grid_Fem_axes(1) % bounds - xBoundsRef)), 0.0_R64, thr=5.1e-5_R64)
  if (allocated(error)) error stop "T_Grid_Fem x element boundaries failure"
  call check(error, size(Grid_Fem_axes(2) % bounds), size(xBoundsRef))
  if (allocated(error)) error stop "T_Grid_Fem y element count failure"
  call check(error, maxval(abs(Grid_Fem_axes(2) % bounds - xBoundsRef)), 0.0_R64, thr=5.1e-5_R64)
  if (allocated(error)) error stop "T_Grid_Fem y element boundaries failure"
  call check(error, size(Grid_Fem_axes(3) % bounds), size(zBoundsRef))
  if (allocated(error)) error stop "T_Grid_Fem z element count failure"
  call check(error, maxval(abs(Grid_Fem_axes(3) % bounds - zBoundsRef)), 0.0_R64, thr=5.1e-5_R64)
  if (allocated(error)) error stop "T_Grid_Fem z element boundaries failure"

  do iCenter = 1, Grid_Fem_nCenters
    do iAxis = 1, 3
      center = Grid_Fem_centers(iAxis, iCenter)
      call check(error, minval(abs(Grid_Fem_axes(iAxis) % bounds - center)) < tiny(1.0_R64), &
                 "center coordinate is not an element boundary")
      if (allocated(error)) error stop "T_Grid_Fem center bound failure"
      call check(error, minval(abs(Grid_Fem_axes(iAxis) % points - center)) < tiny(1.0_R64), &
                 "center coordinate is not a node")
      if (allocated(error)) error stop "T_Grid_Fem center node failure"
    end do
  end do

  !------------------------------------
  ! (b) quadrature of a Gaussian: sum w exp(-r^2) = pi^1.5
  !------------------------------------

  g = Gaussian(1.0_R64, [0.0_R64, 0.0_R64, 0.0_R64], 0.0_R64)

  integral = sum(Grid_Fem_weights * real(g, kind=R64))
  call check(error, integral, PI**1.5_R64, thr=3.5e-8_R64)
  if (allocated(error)) error stop "T_Grid_Fem Gaussian quadrature failure"

  !------------------------------------
  ! (c) kinetic energy of exp(-r^2): -1/2 <g|lap g>/<g|g> = 3/2
  !------------------------------------

  allocate (lapG(Grid_nPoints))
  call Grid_Fem_ApplyLaplacian(lapG, g)

  kinetic = -0.5_R64 * real(Grid_InnerProduct(g, lapG), kind=R64) / real(Grid_InnerProduct(g, g), kind=R64)
  call check(error, kinetic, 1.5_R64, thr=2.0e-6_R64)
  if (allocated(error)) error stop "T_Grid_Fem kinetic energy failure"

  !------------------------------------
  ! (d) self-adjointness in the W metric for two complex off-center functions
  !------------------------------------

  f = Gaussian(0.8_R64, [0.3_R64, -0.2_R64, 0.5_R64], 0.3_R64)
  h = Gaussian(1.1_R64, [-0.4_R64, 0.1_R64, -0.6_R64], -0.5_R64)

  allocate (lapF(Grid_nPoints), lapH(Grid_nPoints))
  call Grid_Fem_ApplyLaplacian(lapF, f)
  call Grid_Fem_ApplyLaplacian(lapH, h)

  lhs = Grid_InnerProduct(f, lapH)
  rhs = Grid_InnerProduct(lapF, h)
  relSymmetry = abs(lhs - rhs) / abs(lhs)
  call check(error, relSymmetry, 0.0_R64, thr=1.5e-13_R64)
  if (allocated(error)) error stop "T_Grid_Fem Laplacian symmetry failure"

  !------------------------------------
  ! (e) eigenbasis: FromEigenbasis(ToEigenbasis(f)) = f and
  ! ToEigenbasis(lap f) = -D ToEigenbasis(f), D(i, j, k) = lx(i) + ly(j) + lz(k)
  !------------------------------------

  allocate (c(Grid_nPoints), back(Grid_nPoints), cLap(Grid_nPoints))
  call Grid_Fem_ToEigenbasis(c, f)
  call Grid_Fem_FromEigenbasis(back, c)

  relRoundTrip = maxval(abs(back - f)) / maxval(abs(f))
  call check(error, relRoundTrip, 0.0_R64, thr=3.0e-14_R64)
  if (allocated(error)) error stop "T_Grid_Fem eigenbasis round trip failure"

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
  call check(error, relDiagonal, 0.0_R64, thr=8.0e-12_R64)
  if (allocated(error)) error stop "T_Grid_Fem eigenbasis diagonalization failure"

  !------------------------------------
  ! (f) shifted solve: (-lap + s) SolveShiftedLaplacian(f, s) = f
  !------------------------------------

  allocate (u(Grid_nPoints), lapU(Grid_nPoints))
  call Grid_Fem_SolveShiftedLaplacian(u, f, shift)
  call Grid_Fem_ApplyLaplacian(lapU, u)

  relResidual = maxval(abs(-lapU + shift * u - f)) / maxval(abs(f))
  call check(error, relResidual, 0.0_R64, thr=2.5e-11_R64)
  if (allocated(error)) error stop "T_Grid_Fem shifted solve failure"

  ! complex shift (as for Crank-Nicolson factors): the real and imaginary
  ! coefficients mix on the diagonal only
  call Grid_Fem_SolveShiftedLaplacian(u, f, complexShift)
  call Grid_Fem_ApplyLaplacian(lapU, u)

  relResidualComplex = maxval(abs(-lapU + complexShift * u - f)) / maxval(abs(f))
  call check(error, relResidualComplex, 0.0_R64, thr=1.5e-12_R64)
  if (allocated(error)) error stop "T_Grid_Fem complex-shift solve failure"

  !------------------------------------
  ! (g) Gram-Schmidt of two overlapping Gaussians with the FEM metric
  !------------------------------------

  allocate (orbSet(Grid_nPoints, 2))
  orbSet(:, 1) = Gaussian(1.0_R64, [0.0_R64, 0.0_R64, 0.7_R64], 0.0_R64)
  orbSet(:, 2) = Gaussian(1.0_R64, [0.0_R64, 0.0_R64, -0.7_R64], 0.0_R64)

  call Grid_Orthonormalize(orbSet)

  ovlp12 = Grid_InnerProduct(orbSet(:, 1), orbSet(:, 2))
  norm11 = Grid_InnerProduct(orbSet(:, 1), orbSet(:, 1))
  norm22 = Grid_InnerProduct(orbSet(:, 2), orbSet(:, 2))

  call check(error, abs(ovlp12), 0.0_R64, thr=5.0e-14_R64)
  if (allocated(error)) error stop "T_Grid_Fem orthonormalize overlap failure"
  call check(error, abs(norm11 - 1.0_R64), 0.0_R64, thr=5.0e-14_R64)
  if (allocated(error)) error stop "T_Grid_Fem orthonormalize norm failure"
  call check(error, abs(norm22 - 1.0_R64), 0.0_R64, thr=5.0e-14_R64)
  if (allocated(error)) error stop "T_Grid_Fem orthonormalize norm failure"

  !------------------------------------
  ! timings (complex fields, averages over nRepeat calls)
  !------------------------------------

  call system_clock(tick0, tickRate)
  do iRepeat = 1, nRepeat
    call Grid_Fem_ApplyLaplacian(lapF, f)
  end do
  call system_clock(tick1)
  tLaplacian = real(tick1 - tick0, R64) / real(tickRate, R64) / nRepeat

  call system_clock(tick0)
  do iRepeat = 1, nRepeat
    call Grid_Fem_ToEigenbasis(c, f)
  end do
  call system_clock(tick1)
  tToEigen = real(tick1 - tick0, R64) / real(tickRate, R64) / nRepeat

  call system_clock(tick0)
  do iRepeat = 1, nRepeat
    call Grid_Fem_SolveShiftedLaplacian(u, f, shift)
  end do
  call system_clock(tick1)
  tShifted = real(tick1 - tick0, R64) / real(tickRate, R64) / nRepeat

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(A, ES12.4)') " (b) sum w exp(-r^2) - pi^1.5                 = ", integral - PI**1.5_R64
  write (*, '(A, F16.10)') " (c) kinetic energy of exp(-r^2)              = ", kinetic
  write (*, '(A, ES12.4)') " (d) |<f|lap h> - <lap f|h>| / |<f|lap h>|     = ", relSymmetry
  write (*, '(A, ES12.4)') " (e) round trip, max rel. error               = ", relRoundTrip
  write (*, '(A, ES12.4)') " (e) ToEigenbasis(lap f) + D c, rel. error    = ", relDiagonal
  write (*, '(A, ES12.4)') " (f) shifted solve, rel. residual             = ", relResidual
  write (*, '(A, ES12.4)') " (f) complex shift 3 - 40i, rel. residual     = ", relResidualComplex
  write (*, '(A, 3ES12.4)') " (g) |<1|2>|, |<1|1> - 1|, |<2|2> - 1|        = ", abs(ovlp12), abs(norm11 - 1.0_R64), &
    abs(norm22 - 1.0_R64)
  write (*, '(A, F9.4, A)') " time ApplyLaplacian (complex)                = ", tLaplacian, " s"
  write (*, '(A, F9.4, A)') " time ToEigenbasis (complex)                  = ", tToEigen, " s"
  write (*, '(A, F9.4, A)') " time SolveShiftedLaplacian (complex)         = ", tShifted, " s"
  print *
  print *, "all fem grid checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> exp(-alpha |r - r0|^2) (1 + i kappa (x - x0)) on the nodes.
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
          field(Grid_Fem_Index(jx, jy, jz)) = exp(-alpha * sum(d**2)) * cmplx(1.0_R64, kappa * d(1), kind=R64)
        end do
      end do
    end do

  end function

end program
