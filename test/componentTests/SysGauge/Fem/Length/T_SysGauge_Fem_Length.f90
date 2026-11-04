!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Component test for the length-gauge laser coupling sysGauge.fem.length,
!> H_E(t) = E(t) (e_pol . r), on a small FEM mesh (box +-5, nLocals 5).
!>
!> Pulse: E0 = 0.05, omega = 0.3, nCycles = 2, cep = 0.4, tStart = 3,
!> polarization (1, 2, -2) -> e_pol = (1, 2, -2)/3.
!>
!> Checks:
!> - flags (time-dependent, body-type-independent), normalized polarization
!> - field amplitude E(t) = E0 sin^2(pi tRel/T) cos(omega (tRel - T/2) + cep)
!>   at times inside the pulse (including the envelope peak) and zero before
!>   and after it
!> - MultiplyWithGaugeOp = E(t) (e_pol . r) orb at every node (oblique
!>   polarization), zero outside the pulse, independent of bt_
!> - Hermiticity in the grid metric: <f|H g> = conjg(<g|H f>)
program T_SysGauge_Fem_Length
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use M_SysGauge
  use M_SysGauge_Fem_Length
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: fieldStrength = 0.05_R64
  real(R64), parameter :: omega = 0.3_R64
  real(R64), parameter :: nCycles = 2.0_R64
  real(R64), parameter :: cep = 0.4_R64
  real(R64), parameter :: tStart = 3.0_R64
  real(R64), parameter :: insideFractions(4) = [0.05_R64, 0.2_R64, 0.5_R64, 0.77_R64]

  integer(I32) :: ix, iy, iz, i, k
  real(R64) :: duration, polRef(3), time, field, fieldRef, errField, errOutside, errPolarization
  real(R64) :: errMultiply, errBodyType, relHermiticity, polDotR, scaleMax
  complex(R64) :: lhs, rhs
  complex(R64), allocatable :: f(:), g(:), hf(:), hg(:), hgBt2(:)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/SysGauge/Fem/Length/T_SysGauge_Fem_Length.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysGauge_Fabricate

  call Grid_Setup
  call SysGauge_Setup

  allocate (hf(Grid_nPoints), hg(Grid_nPoints), hgBt2(Grid_nPoints))

  !==========================================
  call Say_Section("perform action")
  !==========================================

  !------------------------------------
  ! flags and polarization
  !------------------------------------

  call check(error, (.not. SysGauge_timeIndependentQ) .and. SysGauge_bodyTypeIndependentQ, &
             "the laser coupling must be time-dependent and body-type-independent")
  if (allocated(error)) error stop "T_SysGauge_Fem_Length flags failure"

  polRef = [1.0_R64, 2.0_R64, -2.0_R64] / 3.0_R64
  errPolarization = maxval(abs(SysGauge_Fem_Length_polarization - polRef))

  call check(error, errPolarization, 0.0_R64, thr=1.0e-15_R64)
  if (allocated(error)) error stop "T_SysGauge_Fem_Length polarization failure"

  !------------------------------------
  ! field amplitude inside, before and after the pulse
  !------------------------------------

  duration = 2.0_R64 * PI * nCycles / omega

  errField = 0.0_R64
  do k = 1, size(insideFractions)
    time = tStart + insideFractions(k) * duration
    fieldRef = fieldStrength * sin(PI * insideFractions(k))**2 * &
               cos(omega * (insideFractions(k) * duration - 0.5_R64 * duration) + cep)
    errField = max(errField, abs(SysGauge_Fem_Length_FieldAmplitude(time) - fieldRef))
  end do

  call check(error, errField, 0.0_R64, thr=1.0e-15_R64)
  if (allocated(error)) error stop "T_SysGauge_Fem_Length field amplitude failure"

  ! the envelope peak: E0 cos(cep)
  call check(error, SysGauge_Fem_Length_FieldAmplitude(tStart + 0.5_R64 * duration), fieldStrength * cos(cep), &
             thr=1.0e-15_R64)
  if (allocated(error)) error stop "T_SysGauge_Fem_Length field at the envelope peak failure"

  errOutside = max(abs(SysGauge_Fem_Length_FieldAmplitude(0.0_R64)), &
                   abs(SysGauge_Fem_Length_FieldAmplitude(tStart - 0.01_R64)), &
                   abs(SysGauge_Fem_Length_FieldAmplitude(tStart + duration + 0.01_R64)), &
                   abs(SysGauge_Fem_Length_FieldAmplitude(tStart + 3.0_R64 * duration)))

  call check(error, errOutside <= 0.0_R64, "the field must vanish outside the pulse")
  if (allocated(error)) error stop "T_SysGauge_Fem_Length field outside the pulse failure"

  !------------------------------------
  ! MultiplyWithGaugeOp = E(t) (e_pol . r) orb at every node
  !------------------------------------

  f = Gaussian(0.4_R64, [0.3_R64, -0.2_R64, 0.5_R64], 0.3_R64)
  g = Gaussian(0.6_R64, [-0.4_R64, 0.1_R64, -0.6_R64], -0.5_R64)

  time = tStart + 0.3_R64 * duration
  field = SysGauge_Fem_Length_FieldAmplitude(time)

  call SysGauge_MultiplyWithGaugeOp(hf, f, time, bt_=1)

  errMultiply = 0.0_R64
  scaleMax = 0.0_R64
  do iz = 1, Grid_Fem_nz
    do iy = 1, Grid_Fem_ny
      do ix = 1, Grid_Fem_nx
        i = Grid_Fem_Index(ix, iy, iz)
        polDotR = dot_product(polRef, [Grid_Fem_axes(1) % points(ix), Grid_Fem_axes(2) % points(iy), &
                                       Grid_Fem_axes(3) % points(iz)])
        errMultiply = max(errMultiply, abs(hf(i) - field * polDotR * f(i)))
        scaleMax = max(scaleMax, abs(field * polDotR * f(i)))
      end do
    end do
  end do
  errMultiply = errMultiply / scaleMax

  call check(error, errMultiply, 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_SysGauge_Fem_Length MultiplyWithGaugeOp failure"

  ! body type and time outside the pulse
  call SysGauge_MultiplyWithGaugeOp(hg, g, time, bt_=1)
  call SysGauge_MultiplyWithGaugeOp(hgBt2, g, time, bt_=2)
  errBodyType = maxval(abs(hgBt2 - hg))

  call check(error, errBodyType <= 0.0_R64, "the coupling must not depend on the body type")
  if (allocated(error)) error stop "T_SysGauge_Fem_Length body type failure"

  call SysGauge_MultiplyWithGaugeOp(hgBt2, g, tStart + 2.0_R64 * duration)

  call check(error, maxval(abs(hgBt2)) <= 0.0_R64, "the coupling must vanish outside the pulse")
  if (allocated(error)) error stop "T_SysGauge_Fem_Length coupling outside the pulse failure"

  !------------------------------------
  ! Hermiticity in the grid metric
  !------------------------------------

  lhs = Grid_InnerProduct(f, hg)
  rhs = conjg(Grid_InnerProduct(g, hf))
  relHermiticity = abs(lhs - rhs) / abs(lhs)

  call check(error, relHermiticity, 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_SysGauge_Fem_Length Hermiticity failure"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(A, I0)') " Grid_nPoints                                  = ", Grid_nPoints
  write (*, '(A, ES12.4)') " polarization error                            = ", errPolarization
  write (*, '(A, ES12.4)') " field amplitude error inside the pulse        = ", errField
  write (*, '(A, ES12.4)') " field amplitude outside the pulse             = ", errOutside
  write (*, '(A, ES12.4)') " max|H f - E (e_pol . r) f| / max|E (e.r) f|   = ", errMultiply
  write (*, '(A, ES12.4)') " |<f|H g> - conjg<g|H f>| / |<f|H g>|          = ", relHermiticity
  print *
  print *, "all fem length-gauge checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> exp(-alpha |r - r0|^2) (1 + i kappa (x - x0)) on the nodes.
  function Gaussian(alpha, r0, kappa) result(values)
    real(R64), intent(in) :: alpha
    real(R64), intent(in) :: r0(3)
    real(R64), intent(in) :: kappa
    complex(R64), allocatable :: values(:)

    integer(I32) :: jx, jy, jz
    real(R64) :: d(3)

    allocate (values(Grid_nPoints))
    do jz = 1, Grid_Fem_nz
      do jy = 1, Grid_Fem_ny
        do jx = 1, Grid_Fem_nx
          d = [Grid_Fem_axes(1) % points(jx), Grid_Fem_axes(2) % points(jy), Grid_Fem_axes(3) % points(jz)] - r0
          values(Grid_Fem_Index(jx, jy, jz)) = exp(-alpha * sum(d**2)) * cmplx(1.0_R64, kappa * d(1), kind=R64)
        end do
      end do
    end do

  end function

end program
