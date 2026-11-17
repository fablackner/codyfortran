!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Imaginary-time ground state of an orbital-based method (tdhx, mctdhx,
! tdci), shared by the imaginary-time tests of He1d, Ne3d, H2, FermiHubbard
! and BoseHubbard; included after `contains`.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Loads jsonFileName and propagates Method_state with
!> Method_ImagTimeDerivative in output steps of program.outputStep
!> (program.nTimeSteps steps each) until the energy changes by less than
!> program.convThresh in one output step. The orbitals are orthonormalized
!> after every step, the CI vector is normalized as well with
!> normalizeCoeffsQ_ (tdci, mctdhx). Stops unless the final energy agrees
!> with referenceEnergy to threshold. An input without a sysPotential block
!> (lattices) skips SysPotential, one without a configList (coeffs.hubbard)
!> ConfigList. The modules stay set up, so the caller can save the ground
!> state.
subroutine CheckImagTimeGroundState(jsonFileName, referenceEnergy, threshold, normalizeCoeffsQ_)
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
  use M_IntegratorList
  use M_Propagator
  use testdrive, only: check, error_type

  character(len=*), intent(in) :: jsonFileName
  real(R64), intent(in) :: referenceEnergy
  real(R64), intent(in) :: threshold
  logical, intent(in), optional :: normalizeCoeffsQ_

  real(R64) :: energyNew, energyOld, conv, convThresh, outputStep, time, timeStep
  integer(I32) :: iStep, nTimeSteps, innerStep
  logical :: normalizeCoeffsQ, potentialQ, configListQ
  type(error_type), allocatable :: error
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

  normalizeCoeffsQ = .false.
  if (present(normalizeCoeffsQ_)) normalizeCoeffsQ = normalizeCoeffsQ_

  !==========================================
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)
  ! coeffs.hubbard has no configuration list
  configListQ = Json_GetExistence("configList")
  potentialQ = Json_GetExistence("sysPotential")

  call Grid_Fabricate
  call SysKinetic_Fabricate
  if (potentialQ) call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  if (configListQ) call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_ImagTimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Say_Fabricate("program")
  convThresh = Json_Get("program.convThresh", 1e-13_R64)
  outputStep = Json_Get("program.outputStep", 1e-1_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 10)

  call Grid_Setup
  call SysKinetic_Setup
  if (potentialQ) call SysPotential_Setup
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

  iStep = 0
  time = 0.0_R64
  timeStep = outputStep / nTimeSteps
  energyNew = 1e10_R64
  conv = 1e10_R64

  do while (abs(conv) > convThresh)
    iStep = iStep + 1

    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep

      call Orbs_Orthonormalize(Orbs_orbs)
      if (normalizeCoeffsQ) call Coeffs_Normalize(Coeffs_coeffs)
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

  call check(error, energyNew, referenceEnergy, thr=threshold)
  if (allocated(error)) error stop "CheckImagTimeGroundState: energy off the reference"

end subroutine
