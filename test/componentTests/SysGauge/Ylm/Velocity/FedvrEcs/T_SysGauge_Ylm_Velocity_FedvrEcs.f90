!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Component test for the velocity-gauge coupling H_A = A(t) p_z + A(t)²/2
!> on the Ylm FEDVR-ECS grid (complex contour).
!>
!> Checks that the action on a pure (l=0, m=0) channel f(z) = z e^{-λz},
!> analytically continued onto the ECS contour, matches the analytic result:
!> the raising channel (1,0) receives -i A f'(z)/sqrt(3) and the source
!> channel keeps the diagonal A²/2 f(z).
program T_SysGauge_Ylm_Velocity_FedvrEcs
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Ylm
  use M_Grid_Ylm_FedvrEcs
  use M_SysGauge
  use M_SysGauge_Ylm_Velocity
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: i
  real(R64) :: lambda, time, aPot, duration
  real(R64) :: avgError, avgResult
  complex(R64) :: z, analytical
  complex(R64), allocatable :: orb(:), dOrb(:), fLm(:), dOrbLm(:)
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/SysGauge/Ylm/Velocity/FedvrEcs/T_SysGauge_Ylm_Velocity_FedvrEcs.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysGauge_Fabricate

  call Grid_Setup
  call SysGauge_Setup

  allocate (orb(Grid_nPoints), dOrb(Grid_nPoints))
  allocate (fLm(Grid_Ylm_nRadial), dOrbLm(Grid_Ylm_nRadial))

  ! Evaluate the pulse away from the envelope zeros/nodes
  duration = 2.0_R64 * PI * SysGauge_Ylm_Velocity_nCycles / SysGauge_Ylm_Velocity_omega
  time = 0.3_R64 * duration
  aPot = SysGauge_Ylm_Velocity_VectorPotentialAmplitude(time)

  call check(error, abs(aPot) > 1e-3_R64)
  if (allocated(error)) error stop "T_SysGauge_Ylm_Velocity_FedvrEcs failure: vector potential vanishes at test time"

  lambda = 1.0_R64

  ! Initialize test wavefunction f(z) = z e^{-λz} on the (0,0) channel,
  ! analytically continued onto the ECS contour
  orb(:) = 0.0_R64
  do i = 1, Grid_Ylm_nRadial
    z = Grid_Ylm_FedvrEcs_contourPoints(i)
    fLm(i) = z * exp(-lambda * z)
  end do
  call Grid_Ylm_SetLmComponent(orb, 0, 0, fLm)

  call SysGauge_MultiplyWithGaugeOp(dOrb, orb, time, 1)

  ! Raising channel (1,0): -i A c_00 f'(z) with c_00 = 1/sqrt(3),
  ! f'(z) = (1 - λz) e^{-λz}. The included free endpoint at rmax carries the
  ! weak-form boundary term (natural boundary condition), so it is excluded
  ! from the error metric.
  call Grid_Ylm_GetLmComponent(dOrbLm, 1, 0, dOrb)

  avgError = 0.0_R64
  avgResult = 0.0_R64
  do i = 1, Grid_Ylm_nRadial - 1
    z = Grid_Ylm_FedvrEcs_contourPoints(i)
    analytical = -IU * aPot / sqrt(3.0_R64) * (1.0_R64 - lambda * z) * exp(-lambda * z)
    avgError = avgError + abs(dOrbLm(i) - analytical)
    avgResult = avgResult + abs(analytical)
  end do
  write (*, '(A, E12.5)') "Average error (raising channel): ", avgError / avgResult

  call check(error, avgError / avgResult < 1e-5_R64)
  if (allocated(error)) error stop "T_SysGauge_Ylm_Velocity_FedvrEcs failure: raising channel off analytic result"

  ! Source channel (0,0): diagonal A²/2 f(z)
  call Grid_Ylm_GetLmComponent(dOrbLm, 0, 0, dOrb)

  avgError = 0.0_R64
  avgResult = 0.0_R64
  do i = 1, Grid_Ylm_nRadial
    analytical = 0.5_R64 * aPot**2 * fLm(i)
    avgError = avgError + abs(dOrbLm(i) - analytical)
    avgResult = avgResult + abs(analytical)
  end do
  write (*, '(A, E12.5)') "Average error (diagonal A²/2):   ", avgError / avgResult

  call check(error, avgError / avgResult < 1e-12_R64)
  if (allocated(error)) error stop "T_SysGauge_Ylm_Velocity_FedvrEcs failure: diagonal channel off analytic result"

  deallocate (orb, dOrb, fLm, dOrbLm)

  !==========================================
  call Say_Goodbye
  !==========================================

end program
