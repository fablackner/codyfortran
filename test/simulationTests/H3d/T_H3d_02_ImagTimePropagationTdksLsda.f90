!> ctest-requires: libxc
!> Spin-polarized Kohn-Sham ground state of the H atom (one electron, one
!> body type: rho_down = 0, the polarized functional; Slater exchange + VWN5)
!> by imaginary-time propagation of method.mb.orbBased.tdks on a Ylm FEDVR
!> grid, against the NIST LSD reference E = -0.478671 (an independent radial
!> finite-difference solver with the same Libxc functionals gives
!> -0.4786708). The orbital is normalized after every step: the Kohn-Sham
!> potential is nonlinear in the density (normalizing only every output
!> step converges to -0.478138, the fixed point of a drifting norm).
program T_H3d_02_ImagTimePropagationTdksLsda
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_SysXc
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb_OrbBased_Tdks
  use M_IntegratorList
  use M_Propagator
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)
  real(R64) :: energy, energyOld, convThresh, timeStep, outputStep, time
  integer(I32) :: iOut, iStep, nStepsPerOutput, maxOutputs

  call Say_Hello
  call Json_LoadJsonFile(manualFileName_="test/simulationTests/H3d/T_H3d_02_ImagTimePropagationTdksLsda.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call SysXc_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_ImagTimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()
  convThresh = Json_Get("program.convThresh", 1e-11_R64)
  timeStep = Json_Get("program.timeStep", 4e-4_R64)
  outputStep = Json_Get("program.outputStep", 0.5_R64)
  maxOutputs = Json_Get("program.maxOutputs", 200)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call SysXc_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup
  call Orbs_Orthonormalize(Orbs_orbs)

  call check(error, Method_Mb_OrbBased_Tdks_nSpin, 2)
  if (allocated(error)) error stop "T_H3d_02_ImagTimePropagationTdksLsda failure: one body type must be spin-polarized"

  nStepsPerOutput = nint(outputStep / timeStep)
  energy = Method_GetEnergy(0.0_R64)
  time = 0.0_R64
  do iOut = 1, maxOutputs
    ! Normalized after every step: the Kohn-Sham potential is nonlinear in
    ! the density, so the norm must not drift between the steps
    do iStep = 1, nStepsPerOutput
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
      call Orbs_Orthonormalize(Orbs_orbs)
    end do
    energyOld = energy
    energy = Method_GetEnergy(0.0_R64)
    write (*, '(1X, A, F8.3, A, F18.12, A, ES10.2, A, 2F14.10)') "tau ", time, ": E ", energy, ", dE ", &
      energy - energyOld, ", E_H, E_xc ", Method_Mb_OrbBased_Tdks_hartreeEnergy, Method_Mb_OrbBased_Tdks_xcEnergy
    if (abs(energy - energyOld) < convThresh) exit
  end do

  print *
  print '(1X, A, F16.10)', "H LSDA: E = ", energy

  call check(error, abs(energy - energyOld) < convThresh)
  if (allocated(error)) error stop "T_H3d_02_ImagTimePropagationTdksLsda failure: not converged"
  call check(error, energy, -0.478671_R64, thr=2e-6_R64)
  if (allocated(error)) error stop "T_H3d_02_ImagTimePropagationTdksLsda failure: energy (NIST LSD -0.478671)"

  call Say_Goodbye

end program
