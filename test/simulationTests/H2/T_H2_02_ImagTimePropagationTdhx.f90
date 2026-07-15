!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> H2 Hartree-Fock ground state via imaginary time propagation (TDHX) on the
!> prolate-spheroidal grid.
!>
!> Two electrons (spin up/down as two fermionic body types with one orbital
!> each) in the field of two protons at the fixed equilibrium distance
!> R = 1.4, starting from the gerade LCAO 1s sigma_g guess. The converged
!> electronic energy is compared against the Hartree-Fock limit
!>
!>   E_HF(total, R = 1.4) = -1.13362957 (Kolos/Roothaan-quality HF limit)
!>   E_el = E_HF - 1/R    = -1.13362957 - 0.71428571 = -1.84791528
!>
!> (the constant nuclear repulsion 1/R is not part of the electronic
!> Hamiltonian and is excluded from the comparison).
program T_H2_02_ImagTimePropagationTdhx
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Actions
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
  use M_IntegratorList
  use M_Propagator
  use testdrive, only: check, error_type

  implicit none

  real(R64)    :: energyNew = 1e10_R64
  real(R64)    :: conv = 1e10_R64
  real(R64)    :: energyOld
  real(R64)    :: convThresh
  real(R64)    :: outputStep
  real(R64)    :: time
  real(R64)    :: timeStep
  integer(I32) :: iStep, nTimeSteps, innerStep
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/simulationTests/H2/T_H2_02_ImagTimePropagationTdhx.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  IntegratorListInput(1) % TimeDerivative => Actions_ImagTimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Say_Fabricate("program")
  convThresh = Json_Get("program.convThresh", 1e-13_R64)
  outputStep = Json_Get("program.outputStep", 1e-1_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 10)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  iStep = 0
  time = 0.0_R64
  timeStep = outputStep / nTimeSteps

  do while (abs(conv) > convThresh)
    iStep = iStep + 1

    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep

      call Orbs_Orthonormalize(Orbs_orbs)
    end do

    energyOld = energyNew
    energyNew = Method_GetEnergy(time)

    conv = energyNew - energyOld
    write (*, *) "convergence: ", energyNew, conv, " time:", time
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "final energy: ", energyNew

  ! Hartree-Fock limit of the electronic energy at R = 1.4 (see header);
  ! the converged grid value -1.84791526025 agrees with it to ~2e-8
  call check(error, energyNew, -1.8479152602_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_H2_02_ImagTimePropagationTdhx failure"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
