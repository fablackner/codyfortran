! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Spin-restricted TDHX imaginary-time propagation for 1D helium.
!>
!> Identical physics to T_He1d_03_ImagTimePropagationTdhx, but with
!> orbs.restrictedQ = true: only the shared spatial orbital set is stored and
!> propagated, and the mean-field (Hartree) term uses the spin-summed density.
!> The converged Hartree-Fock energy must therefore match the unrestricted
!> reference value.
program T_He1d_08_ImagTimePropagationTdhxRestricted
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

  jsonFileName = "test/simulationTests/He1d/T_He1d_08_ImagTimePropagationTdhxRestricted.json"
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
  time = 0.0_R64  ! Initialize simulation time
  timeStep = outputStep / nTimeSteps  ! Calculate step size for each iteration

  do while (abs(conv) > convThresh)
    iStep = iStep + 1

    ! Explicit loop to execute propagation nTimeSteps times
    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep  ! Update current time

      call Orbs_Orthonormalize(Orbs_orbs)
    end do

    energyOld = energyNew
    energyNew = Method_GetEnergy(time)  ! Use current time for energy calculation

    conv = energyNew - energyOld
    write (*, *) "convergence: ", energyNew, conv, " time:", time
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "final energy: ", energyNew

  ! Same reference energy as the unrestricted T_He1d_03 test
  call check(error, energyNew, -2.4819265615333204_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_He1d_ImagTimePropagationTdhxRestricted failure"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
