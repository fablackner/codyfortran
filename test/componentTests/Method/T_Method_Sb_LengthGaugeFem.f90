!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> method.sb with the length-gauge laser of grid.fem (sysGauge.fem.length):
!> the single-body time derivative contains the laser coupling. H atom on a
!> tiny mesh (box ±6, nLocals 4, hMin 0.5), 1s Slater start (orbsInit.fem.sto).
!>
!> Pulse along z, E(t) = F0 sin^2(pi tRel/T) cos(omega (tRel - T/2) + cep),
!> tRel = t - tStart, T = 2 pi nCycles/omega, zero outside [tStart, tStart + T];
!> at the envelope peak t = tStart + T/2 (cep = 0) E = F0 = 0.05.
!>
!> Checks:
!> - dpsi/dt at the envelope peak = -i (T + V + E(t) z) psi, with E(t) from the
!>   explicit formula above and T, V, z applied here (relative 1e-13)
!> - the laser term is a sizable part of H psi there (the check above is not
!>   trivially met)
!> - before the pulse (t = 0 < tStart) dpsi/dt = -i (T + V) psi
program T_Method_Sb_LengthGaugeFem
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use M_SysKinetic
  use M_SysPotential
  use M_SysGauge
  use M_OrbsInit
  use M_Method
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  complex(R64), allocatable :: psi(:), dPsi(:), hPsi(:), zPsi(:), potential(:)
  real(R64) :: fieldStrength, omega, nCycles, cep, tStart, duration, tPeak, fieldPeak
  real(R64) :: devPeak, devBefore, fieldShare
  integer(I32) :: ix, iy, iz, node

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Method/T_Method_Sb_LengthGaugeFem.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysGauge_Fabricate
  call OrbsInit_Fabricate
  call Method_Fabricate

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysGauge_Setup
  call OrbsInit_Setup
  call Method_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  ! Envelope peak and the field there, from the pulse parameters
  fieldStrength = Json_Get("sysGauge.fem.length.fieldStrength", 0.0_R64)
  omega = Json_Get("sysGauge.fem.length.omega", 1.0_R64)
  nCycles = Json_Get("sysGauge.fem.length.nCycles", 1.0_R64)
  cep = Json_Get("sysGauge.fem.length.cep", 0.0_R64)
  tStart = Json_Get("sysGauge.fem.length.tStart", 0.0_R64)
  duration = 2.0_R64 * PI * nCycles / omega
  tPeak = tStart + 0.5_R64 * duration
  fieldPeak = fieldStrength * sin(PI * (tPeak - tStart) / duration)**2 * &
              cos(omega * (tPeak - tStart - 0.5_R64 * duration) + cep)

  ! (T + V) psi and z psi of the 1s start, applied here
  psi = Method_state
  allocate (dPsi, hPsi, zPsi, mold=psi)
  call SysKinetic_MultiplyWithKineticOp(hPsi, psi, tPeak)
  call SysPotential_FillExternalPotential(potential, tPeak)
  hPsi(:) = hPsi(:) + potential(:) * psi(:)
  do iz = 1, Grid_Fem_nz
    do iy = 1, Grid_Fem_ny
      do ix = 1, Grid_Fem_nx
        node = Grid_Fem_Index(ix, iy, iz)
        zPsi(node) = Grid_Fem_axes(3) % points(iz) * psi(node)
      end do
    end do
  end do

  ! At the envelope peak
  call Method_TimeDerivative(dPsi, Method_state, tPeak)
  devPeak = MetricNorm(dPsi + IU * (hPsi + fieldPeak * zPsi)) / MetricNorm(hPsi + fieldPeak * zPsi)
  fieldShare = MetricNorm(fieldPeak * zPsi) / MetricNorm(hPsi + fieldPeak * zPsi)

  ! Before the pulse
  call Method_TimeDerivative(dPsi, Method_state, 0.0_R64)
  devBefore = MetricNorm(dPsi + IU * hPsi) / MetricNorm(hPsi)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(1X, A, I0)') "mesh N = ", Grid_nPoints
  print *, "envelope peak t, E(t):                        ", tPeak, fieldPeak
  print *, "|dpsi/dt + i (T + V + E z) psi| / |H psi|:    ", devPeak
  print *, "|E z psi| / |H psi| at the peak:              ", fieldShare
  print *, "|dpsi/dt + i (T + V) psi| / |H psi|, t = 0:   ", devBefore

  call check(error, abs(fieldPeak - fieldStrength) < 1e-14_R64)
  if (allocated(error)) error stop "T_Method_Sb_LengthGaugeFem failure: field at the envelope peak is not F0"

  call check(error, devPeak < 1e-13_R64)
  if (allocated(error)) error stop "T_Method_Sb_LengthGaugeFem failure: dpsi/dt lacks -i E(t) z psi at the peak"

  call check(error, fieldShare > 1e-2_R64)
  if (allocated(error)) error stop "T_Method_Sb_LengthGaugeFem failure: laser term too small to be tested"

  call check(error, devBefore < 1e-13_R64)
  if (allocated(error)) error stop "T_Method_Sb_LengthGaugeFem failure: dpsi/dt before the pulse is not -i (T + V) psi"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Norm in the grid metric, sqrt(<f|f>).
  function MetricNorm(f) result(res)
    complex(R64), intent(in), contiguous :: f(:)
    real(R64) :: res

    res = sqrt(real(Grid_InnerProduct(f, f), kind=R64))

  end function

end program
