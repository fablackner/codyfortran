!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Real-time TDCI on an open chain of 5 sites (2 + 2 fermions, hopping 1,
! U = 2, on-site energies 0.3 -0.1 0 0.2 -0.4) from the configuration with
! both spins on sites 1 and 2, which is symmetric under the exchange of the
! spins. Shared by T_FermiHubbard_10 (coeffs.generic), T_FermiHubbard_11
! (coeffs.hubbard.noSpinSym) and T_FermiHubbard_12 (coeffs.hubbard.plusSpinSym);
! included after `contains`. The reference values come from coeffs.generic.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Loads jsonFileName, propagates the unary state with RK4 steps of
!> program.timeStep to program.endTime and stops unless the energy is
!> conserved and the site occupations <n_up + n_dn> and double occupations
!> <n_up n_dn> at the end time agree with the coeffs.generic reference.
subroutine CheckRealTimeTdci(jsonFileName)
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator
  use testdrive, only: check, error_type

  character(len=*), intent(in) :: jsonFileName

  ! coeffs.generic at t = 2 (from 2, 2, 0, 0, 0 and 1, 1, 0, 0, 0 at t = 0)
  real(R64), parameter :: occupationRef(5) = [1.2534476939255705_R64, 1.1191016437652728_R64, 7.6550057985958819e-1_R64, &
                                              4.9106164643206307e-1_R64, 3.7088825642284934e-1_R64]
  real(R64), parameter :: doubleRef(5) = [5.7749218574667460e-1_R64, 4.5958576361803372e-1_R64, 2.6966245218219964e-1_R64, &
                                          1.6268046620067075e-1_R64, 8.6580200675616892e-2_R64]

  real(R64) :: time, timeStep, endTime, energy0, energy
  real(R64) :: occupation(5), double(5)
  integer(I32) :: iStep, nSteps, s, oUp, oDn
  logical :: configListQ
  type(error_type), allocatable :: error
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

  !==========================================
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)
  ! coeffs.hubbard has no configuration list
  configListQ = Json_GetExistence("configList")

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  if (configListQ) call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_TimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Say_Fabricate("program")
  endTime = Json_Get("program.endTime", 1.0_R64)
  timeStep = Json_Get("program.timeStep", 0.01_R64)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  if (configListQ) call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  energy0 = Method_GetEnergy(0.0_R64)
  nSteps = nint(endTime / timeStep)
  time = 0.0_R64
  do iStep = 1, nSteps
    call Propagator_Propagate(Method_state, time, time + timeStep)
    time = time + timeStep
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  ! the energy fills Method_Mb_OrbBased_rdm1 and _rdm2 of the final state
  energy = Method_GetEnergy(time)
  oUp = Method_Mb_OrbBased_nOrbsStart(1) - 1
  oDn = Method_Mb_OrbBased_nOrbsStart(2) - 1
  do s = 1, 5
    occupation(s) = real(Method_Mb_OrbBased_rdm1(oUp + s, oUp + s) + Method_Mb_OrbBased_rdm1(oDn + s, oDn + s), kind=R64)
    double(s) = real(Method_Mb_OrbBased_rdm2(oUp + s, oDn + s, oUp + s, oDn + s), kind=R64)
  end do

  print *
  print *, "energy start, end:     ", energy0, energy
  write (*, '(A,5ES24.16)') " site occupations:     ", occupation
  write (*, '(A,5ES24.16)') " double occupations:   ", double

  ! RK4 with dt = 0.01 changes the energy 4.4 by 2.5e-7 over 2 a.u.
  call check(error, energy, energy0, thr=1e-6_R64)
  if (allocated(error)) error stop "CheckRealTimeTdci: energy not conserved"
  call check(error, maxval(abs(occupation - occupationRef)) < 1e-10_R64)
  if (allocated(error)) error stop "CheckRealTimeTdci: site occupations differ from the coeffs.generic reference"
  call check(error, maxval(abs(double - doubleRef)) < 1e-10_R64)
  if (allocated(error)) error stop "CheckRealTimeTdci: double occupations differ from the coeffs.generic reference"

end subroutine
