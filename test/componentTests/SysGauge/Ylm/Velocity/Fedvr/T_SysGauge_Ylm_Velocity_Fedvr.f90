!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Component test for the velocity-gauge coupling H_A = A(t) p_z + A(t)²/2
!> on the Ylm FEDVR grid.
!>
!> Checks:
!> - the action on a pure (l=0, m=0) channel f(r) = r e^{-λr} matches the
!>   analytic result: the raising channel (1,0) receives -i A f'(r)/sqrt(3)
!>   and the source channel keeps the diagonal A²/2 f(r)
!> - channels not coupled by a z-polarized field stay empty
!> - the operator is Hermitian in the FEDVR-weighted metric
program T_SysGauge_Ylm_Velocity_Fedvr
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Ylm
  use M_SysGauge
  use M_SysGauge_Ylm_Velocity
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: i
  real(R64) :: r, lambda, time, aPot, duration
  real(R64) :: avgError, avgResult
  complex(R64), allocatable :: orb(:), dOrb(:), orbBra(:), dOrbBra(:)
  complex(R64), allocatable :: fLm(:), dOrbLm(:)
  complex(R64) :: analytical, braKet, ketBra
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/SysGauge/Ylm/Velocity/Fedvr/T_SysGauge_Ylm_Velocity_Fedvr.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysGauge_Fabricate

  call Grid_Setup
  call SysGauge_Setup

  allocate (orb(Grid_nPoints), dOrb(Grid_nPoints))
  allocate (orbBra(Grid_nPoints), dOrbBra(Grid_nPoints))
  allocate (fLm(Grid_Ylm_nRadial), dOrbLm(Grid_Ylm_nRadial))

  ! Evaluate the pulse away from the envelope zeros/nodes
  duration = 2.0_R64 * PI * SysGauge_Ylm_Velocity_nCycles / SysGauge_Ylm_Velocity_omega
  time = 0.3_R64 * duration
  aPot = SysGauge_Ylm_Velocity_VectorPotentialAmplitude(time)

  call check(error, abs(aPot) > 1e-3_R64)
  if (allocated(error)) error stop "T_SysGauge_Ylm_Velocity_Fedvr failure: vector potential vanishes at test time"

  !==========================================
  call Say_Section("analytic action on a (0,0) channel")
  !==========================================

  lambda = 2.0_R64

  ! Initialize test wavefunction f(r) = r e^{-λr} on the (0,0) channel
  orb(:) = 0.0_R64
  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    fLm(i) = r * exp(-lambda * r)
  end do
  call Grid_Ylm_SetLmComponent(orb, 0, 0, fLm)

  call SysGauge_MultiplyWithGaugeOp(dOrb, orb, time, 1)

  ! Raising channel (1,0): -i A c_00 f'(r) with c_00 = 1/sqrt(3),
  ! f'(r) = (1 - λr) e^{-λr}
  call Grid_Ylm_GetLmComponent(dOrbLm, 1, 0, dOrb)

  avgError = 0.0_R64
  avgResult = 0.0_R64
  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    analytical = -IU * aPot / sqrt(3.0_R64) * (1.0_R64 - lambda * r) * exp(-lambda * r)
    avgError = avgError + abs(dOrbLm(i) - analytical)
    avgResult = avgResult + abs(analytical)
  end do
  write (*, '(A, E12.5)') "Average error (raising channel): ", avgError / avgResult

  call check(error, avgError / avgResult < 1e-5_R64)
  if (allocated(error)) error stop "T_SysGauge_Ylm_Velocity_Fedvr failure: raising channel off analytic result"

  ! Source channel (0,0): diagonal A²/2 f(r)
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
  if (allocated(error)) error stop "T_SysGauge_Ylm_Velocity_Fedvr failure: diagonal channel off analytic result"

  ! A z-polarized field couples only Δl = ±1, Δm = 0: (2,0) must stay empty
  call Grid_Ylm_GetLmComponent(dOrbLm, 2, 0, dOrb)

  call check(error, all(abs(dOrbLm) < 1e-14_R64))
  if (allocated(error)) error stop "T_SysGauge_Ylm_Velocity_Fedvr failure: uncoupled channel is populated"

  !==========================================
  call Say_Section("hermiticity in the weighted metric")
  !==========================================

  ! Bra and ket with several coupled channels; the decay is chosen fast
  ! enough that the r_max boundary term of the first derivative (the grid
  ! includes the r_max endpoint) stays below the check tolerance
  orbBra(:) = 0.0_R64
  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    fLm(i) = r * exp(-2.0_R64 * r)
  end do
  call Grid_Ylm_SetLmComponent(orbBra, 0, 0, fLm)
  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    fLm(i) = r**2 * exp(-2.3_R64 * r)
  end do
  call Grid_Ylm_SetLmComponent(orbBra, 1, 0, fLm)

  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    fLm(i) = r * exp(-1.8_R64 * r)
  end do
  call Grid_Ylm_SetLmComponent(orb, 0, 0, fLm)
  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    fLm(i) = r**2 * exp(-2.1_R64 * r)
  end do
  call Grid_Ylm_SetLmComponent(orb, 2, 0, fLm)

  call SysGauge_MultiplyWithGaugeOp(dOrb, orb, time, 1)
  call SysGauge_MultiplyWithGaugeOp(dOrbBra, orbBra, time, 1)

  braKet = Grid_InnerProduct(orbBra, dOrb)
  ketBra = Grid_InnerProduct(orb, dOrbBra)

  write (*, '(A, 2E12.5)') "Hermiticity defect:              ", abs(braKet - conjg(ketBra))

  call check(error, abs(braKet - conjg(ketBra)) < 1e-12_R64)
  if (allocated(error)) error stop "T_SysGauge_Ylm_Velocity_Fedvr failure: operator not Hermitian in weighted metric"

  deallocate (orb, dOrb, orbBra, dOrbBra, fLm, dOrbLm)

  !==========================================
  call Say_Goodbye
  !==========================================

end program
