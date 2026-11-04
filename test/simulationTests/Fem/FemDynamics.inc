!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Helpers shared by the real-time tests T_Fem_05 .. T_Fem_07, T_Fem_09 and
! T_Fem_10 (included after `contains`; every procedure here is used by each of
! them).
!
! - SetupFemDynamics:  fabricate and set up the run of a JSON file (restricted
!                      MCTDHX with the orbital of T_Fem_03 loaded, split-step
!                      propagator, exact field-free part by the Cayley step)
! - FemObservables:    energy, particle number and dipole of Method_state
! - ElapsedSeconds:    wall-clock time
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fabricates and sets up the run of jsonFileName in the order of the other
!> real-time tests (Ne3dMctdhxLaser.inc): split-step propagator with the
!> linear part (gauge coupling only, h0 exact) in integrator 1 and the
!> interaction part in integrator 2. The loaded orbitals are orthonormalized
!> after Method_Setup (OrbsInit never does).
subroutine SetupFemDynamics(jsonFileName)
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysGauge
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

  character(len=*), intent(in) :: jsonFileName

  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)
  logical :: gaugeQ

  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)
  gaugeQ = Json_GetExistence("sysGauge")

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  if (gaugeQ) call SysGauge_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  if (gaugeQ) call SysGauge_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup

  call Orbs_Orthonormalize(Orbs_orbs)

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Electronic energy (Method_GetEnergy, without the nuclear repulsion),
!> particle number N = Tr(rho S) and dipole (<x>, <y>, <z>) of Method_state;
!> the dipole is appended to dipoleFile ("time <x> <y> <z>") if it is not empty.
subroutine FemObservables(energy, norm, dipole, time, dipoleFile)
  use M_Utils_PrinterObservableFem
  use M_Orbs
  use M_Method
  use M_Method_Mb_OrbBased

  real(R64), intent(out) :: energy
  real(R64), intent(out) :: norm
  real(R64), intent(out) :: dipole(3)
  real(R64), intent(in) :: time
  character(len=*), intent(in) :: dipoleFile

  complex(R64), allocatable :: rdm1(:, :)

  energy = Method_GetEnergy(time)
  call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
  call PrinterObservableFem_DumpDipole(rdm1, Orbs_orbs, dipoleFile, .false., time, dipole_=dipole)
  call PrinterObservableFem_DumpNorm(rdm1, Orbs_orbs, "", .false., time, norm_=norm)

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Wall-clock seconds since count0 (system_clock counts).
function ElapsedSeconds(count0) result(res)
  integer(I64), intent(in) :: count0
  real(R64) :: res

  integer(I64) :: count1, countRate

  call system_clock(count1, countRate)
  res = real(count1 - count0, kind=R64) / real(countRate, kind=R64)

end function
