!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Component test for the spectral-element kinetic operator T = -lap/(2m)
!> (sysKinetic.fem.laplacian) and its exact shifted inverse
!> (SysKinetic_SolveShifted) on a mesh graded toward the origin: box +-8,
!> nLocals 6, one center at the origin, hMin 0.2, growth 1.6, hMax 2
!> (n = 69 per axis, N = 328509); bodyMass = [1, 2.5].
!>
!> Checks:
!> - flags: time-independent, body-type-dependent (different masses), the
!>   optional SolveShifted hook is bound
!> - Hermiticity in the grid metric: <f|T g> = conjg(<g|T f>)
!> - <g|T|g>/<g|g> = 3a/2 for g = exp(-a r^2), a = 1 and 0.5
!> - mass scaling: T(bt = 2) = T(bt = 1)/2.5, absent bt_ = body type 1
!> - SolveShifted inverts T + s for the real shift s = 0.7 and the complex
!>   Crank-Nicolson shift s = -2i/dt (dt = 0.05), checked by applying T
!>   explicitly: max|(T + s) x - b| / max|b| ~ 1e-10; a complex shift with
!>   zero imaginary part equals the real-shift solve
!> - bodyMass with a single entry (refabricated with [1.7]): body-type
!>   independent, the mass of every body type (bt_ = 3, beyond the array,
!>   and an absent bt_ give T(bt = 1); <g|T|g>/<g|g> = 1.5/1.7; SolveShifted
!>   for bt_ = 3 inverts T + s); several equal masses [1, 1]: body-type
!>   independent (a body type beyond several masses stops:
!>   T_SysKinetic_Fem_Laplacian_MissingMass)
program T_SysKinetic_Fem_Laplacian
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use M_SysKinetic
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: massRatio = 2.5_R64
  real(R64), parameter :: singleMass = 1.7_R64
  real(R64), parameter :: realShift = 0.7_R64
  real(R64), parameter :: dt = 0.05_R64
  integer(I32), parameter :: nRepeat = 5

  integer(I32) :: iRepeat
  integer(I64) :: tick0, tick1, tickRate
  real(R64) :: relHermiticity, kinetic1, kinetic05, kineticMass2, relMassScaling
  real(R64) :: resRealShift, resComplexShift, resComplexShiftMass2, diffZeroImag
  real(R64) :: relSingleMass, kineticSingle, resSingleMass
  real(R64) :: tKinetic, tSolveReal, tSolveComplex
  complex(R64) :: lhs, rhs, cnShift
  complex(R64), allocatable :: f(:), g(:), tf(:), tg(:), tg2(:), tgDefault(:), x(:), x2(:), tx(:)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/SysKinetic/Fem/Laplacian/T_SysKinetic_Fem_Laplacian.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate

  call Grid_Setup
  call SysKinetic_Setup

  allocate (tf(Grid_nPoints), tg(Grid_nPoints), tg2(Grid_nPoints), tgDefault(Grid_nPoints))
  allocate (x(Grid_nPoints), x2(Grid_nPoints), tx(Grid_nPoints))

  !==========================================
  call Say_Section("perform action")
  !==========================================

  !------------------------------------
  ! flags and hooks
  !------------------------------------

  call check(error, SysKinetic_timeIndependentQ, "kinetic operator must be time-independent")
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian timeIndependentQ failure"
  call check(error, .not. SysKinetic_bodyTypeIndependentQ, "different masses must make T body-type-dependent")
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian bodyTypeIndependentQ failure"
  call check(error, associated(SysKinetic_SolveShifted), "SysKinetic_SolveShifted must be bound")
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian SolveShifted binding failure"

  !------------------------------------
  ! Hermiticity in the grid metric for two complex off-center functions
  !------------------------------------

  f = Gaussian(0.8_R64, [0.3_R64, -0.2_R64, 0.5_R64], 0.3_R64)
  g = Gaussian(1.1_R64, [-0.4_R64, 0.1_R64, -0.6_R64], -0.5_R64)

  call SysKinetic_MultiplyWithKineticOp(tf, f, 0.0_R64, bt_=1)
  call SysKinetic_MultiplyWithKineticOp(tg, g, 0.0_R64, bt_=1)

  lhs = Grid_InnerProduct(f, tg)
  rhs = conjg(Grid_InnerProduct(g, tf))
  relHermiticity = abs(lhs - rhs) / abs(lhs)

  call check(error, relHermiticity, 0.0_R64, thr=1.0e-12_R64)
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian Hermiticity failure"

  !------------------------------------
  ! kinetic energy of exp(-a r^2): 3a/2 (mass 1)
  !------------------------------------

  g = Gaussian(1.0_R64, [0.0_R64, 0.0_R64, 0.0_R64], 0.0_R64)
  call SysKinetic_MultiplyWithKineticOp(tg, g, 0.0_R64, bt_=1)
  kinetic1 = real(Grid_InnerProduct(g, tg), kind=R64) / real(Grid_InnerProduct(g, g), kind=R64)

  call check(error, kinetic1, 1.5_R64, thr=2.0e-6_R64)
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian kinetic energy (a = 1) failure"

  f = Gaussian(0.5_R64, [0.0_R64, 0.0_R64, 0.0_R64], 0.0_R64)
  call SysKinetic_MultiplyWithKineticOp(tf, f, 0.0_R64, bt_=1)
  kinetic05 = real(Grid_InnerProduct(f, tf), kind=R64) / real(Grid_InnerProduct(f, f), kind=R64)

  call check(error, kinetic05, 0.75_R64, thr=3.0e-7_R64)
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian kinetic energy (a = 0.5) failure"

  !------------------------------------
  ! mass scaling: T(bt = 2) = T(bt = 1) / 2.5; absent bt_ = body type 1
  !------------------------------------

  call SysKinetic_MultiplyWithKineticOp(tg2, g, 0.0_R64, bt_=2)
  call SysKinetic_MultiplyWithKineticOp(tgDefault, g, 0.0_R64)

  relMassScaling = maxval(abs(massRatio * tg2 - tg)) / maxval(abs(tg))
  kineticMass2 = real(Grid_InnerProduct(g, tg2), kind=R64) / real(Grid_InnerProduct(g, g), kind=R64)

  call check(error, relMassScaling, 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian mass scaling failure"
  call check(error, kineticMass2, 1.5_R64 / massRatio, thr=1.0e-6_R64)
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian kinetic energy (mass 2.5) failure"
  call check(error, maxval(abs(tgDefault - tg)) <= 0.0_R64, "absent bt_ must give the mass of body type 1")
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian default body type failure"

  !------------------------------------
  ! exact shifted solves: (T + s) x = b, checked with T applied explicitly
  !------------------------------------

  f = Gaussian(0.8_R64, [0.3_R64, -0.2_R64, 0.5_R64], 0.3_R64)

  ! real shift (preconditioner use)
  call SysKinetic_SolveShifted(x, f, cmplx(realShift, 0.0_R64, kind=R64), bt_=1)
  call SysKinetic_MultiplyWithKineticOp(tx, x, 0.0_R64, bt_=1)
  resRealShift = maxval(abs(tx + realShift * x - f)) / maxval(abs(f))

  call check(error, resRealShift, 0.0_R64, thr=1.0e-10_R64)
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian real-shift solve failure"

  ! the complex interface with zero imaginary part takes the real path
  call Grid_Fem_SolveShiftedLaplacian(x2, f, 2.0_R64 * realShift)
  x2 = 2.0_R64 * x2
  diffZeroImag = maxval(abs(x2 - x))

  call check(error, diffZeroImag, 0.0_R64, thr=1.0e-15_R64 * maxval(abs(x)))
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian zero-imaginary shift failure"

  ! complex Crank-Nicolson shift 2/(i dt) = -2i/dt
  cnShift = cmplx(0.0_R64, -2.0_R64 / dt, kind=R64)

  call SysKinetic_SolveShifted(x, f, cnShift, bt_=1)
  call SysKinetic_MultiplyWithKineticOp(tx, x, 0.0_R64, bt_=1)
  resComplexShift = maxval(abs(tx + cnShift * x - f)) / maxval(abs(f))

  call check(error, resComplexShift, 0.0_R64, thr=1.0e-10_R64)
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian complex-shift solve failure"

  ! the same for body type 2 (mass 2.5) and a real right-hand side
  g = Gaussian(1.0_R64, [0.2_R64, 0.0_R64, -0.3_R64], 0.0_R64)

  call SysKinetic_SolveShifted(x, g, cnShift, bt_=2)
  call SysKinetic_MultiplyWithKineticOp(tx, x, 0.0_R64, bt_=2)
  resComplexShiftMass2 = maxval(abs(tx + cnShift * x - g)) / maxval(abs(g))

  call check(error, resComplexShiftMass2, 0.0_R64, thr=1.0e-10_R64)
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian complex-shift solve (mass 2.5) failure"

  !------------------------------------
  ! one bodyMass entry: the mass of every body type
  !------------------------------------

  call jsonFile % add("sysKinetic.fem.laplacian.bodyMass", [singleMass])
  call SysKinetic_Fabricate

  call check(error, SysKinetic_bodyTypeIndependentQ, "a single mass must make T body-type-independent")
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian single-mass flag failure"

  g = Gaussian(1.0_R64, [0.0_R64, 0.0_R64, 0.0_R64], 0.0_R64)
  call SysKinetic_MultiplyWithKineticOp(tg, g, 0.0_R64, bt_=1)
  call SysKinetic_MultiplyWithKineticOp(tg2, g, 0.0_R64, bt_=3)
  call SysKinetic_MultiplyWithKineticOp(tgDefault, g, 0.0_R64)
  relSingleMass = max(maxval(abs(tg2 - tg)), maxval(abs(tgDefault - tg))) / maxval(abs(tg))
  kineticSingle = real(Grid_InnerProduct(g, tg), kind=R64) / real(Grid_InnerProduct(g, g), kind=R64)

  call check(error, relSingleMass <= 0.0_R64, "every body type must get the single mass")
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian single mass for every body type failure"
  call check(error, kineticSingle, 1.5_R64 / singleMass, thr=2.0e-6_R64)
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian kinetic energy (single mass 1.7) failure"

  call SysKinetic_SolveShifted(x, f, cmplx(realShift, 0.0_R64, kind=R64), bt_=3)
  call SysKinetic_MultiplyWithKineticOp(tx, x, 0.0_R64, bt_=3)
  resSingleMass = maxval(abs(tx + realShift * x - f)) / maxval(abs(f))

  call check(error, resSingleMass, 0.0_R64, thr=1.0e-10_R64)
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian real-shift solve (single mass, bt 3) failure"

  ! several equal masses: body-type independent
  call jsonFile % add("sysKinetic.fem.laplacian.bodyMass", [1.0_R64, 1.0_R64])
  call SysKinetic_Fabricate

  call check(error, SysKinetic_bodyTypeIndependentQ, "equal masses must make T body-type-independent")
  if (allocated(error)) error stop "T_SysKinetic_Fem_Laplacian equal-masses flag failure"

  !------------------------------------
  ! timings (complex fields, averages over nRepeat calls)
  !------------------------------------

  call system_clock(tick0, tickRate)
  do iRepeat = 1, nRepeat
    call SysKinetic_MultiplyWithKineticOp(tf, f, 0.0_R64, bt_=1)
  end do
  call system_clock(tick1)
  tKinetic = real(tick1 - tick0, R64) / real(tickRate, R64) / nRepeat

  call system_clock(tick0)
  do iRepeat = 1, nRepeat
    call SysKinetic_SolveShifted(x, f, cmplx(realShift, 0.0_R64, kind=R64), bt_=1)
  end do
  call system_clock(tick1)
  tSolveReal = real(tick1 - tick0, R64) / real(tickRate, R64) / nRepeat

  call system_clock(tick0)
  do iRepeat = 1, nRepeat
    call SysKinetic_SolveShifted(x, f, cnShift, bt_=1)
  end do
  call system_clock(tick1)
  tSolveComplex = real(tick1 - tick0, R64) / real(tickRate, R64) / nRepeat

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(A, I0)') " Grid_nPoints                                   = ", Grid_nPoints
  write (*, '(A, ES12.4)') " |<f|Tg> - conjg<g|Tf>| / |<f|Tg>|             = ", relHermiticity
  write (*, '(A, F16.10)') " <g|T|g>/<g|g>, a = 1   (exact 1.5)            = ", kinetic1
  write (*, '(A, F16.10)') " <g|T|g>/<g|g>, a = 0.5 (exact 0.75)           = ", kinetic05
  write (*, '(A, F16.10)') " <g|T|g>/<g|g>, a = 1, mass 2.5 (exact 0.6)    = ", kineticMass2
  write (*, '(A, ES12.4)') " max|2.5 T(bt=2) g - T(bt=1) g| / max|T g|     = ", relMassScaling
  write (*, '(A, ES12.4)') " SolveShifted s = 0.7, rel. residual           = ", resRealShift
  write (*, '(A, ES12.4)') " complex path, zero imaginary shift, max diff  = ", diffZeroImag
  write (*, '(A, ES12.4)') " SolveShifted s = -2i/dt, rel. residual        = ", resComplexShift
  write (*, '(A, ES12.4)') " SolveShifted s = -2i/dt, mass 2.5, residual   = ", resComplexShiftMass2
  write (*, '(A, ES12.4)') " single mass 1.7: max|T(bt=3) g - T(bt=1) g|    = ", relSingleMass
  write (*, '(A, F16.10)') " single mass 1.7: <g|T|g>/<g|g> (exact 0.88235) = ", kineticSingle
  write (*, '(A, ES12.4)') " single mass 1.7: SolveShifted bt = 3, residual = ", resSingleMass
  write (*, '(A, F9.4, A)') " time MultiplyWithKineticOp (complex)          = ", tKinetic, " s"
  write (*, '(A, F9.4, A)') " time SolveShifted, real shift (complex f)     = ", tSolveReal, " s"
  write (*, '(A, F9.4, A)') " time SolveShifted, complex shift (complex f)  = ", tSolveComplex, " s"
  print *
  print *, "all fem kinetic checks passed"

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
