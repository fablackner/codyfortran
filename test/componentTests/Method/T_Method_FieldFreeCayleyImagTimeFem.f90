!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Expected error: imaginary time with the Cayley exact part of grid.fem.
!> The Crank-Nicolson step does not damp the high modes, so
!> Method_Mb_OrbBased_PropagateFieldFreeCayley refuses imaginary time and
!> points to the SCF ground state. H atom on the mesh of
!> T_Method_Sb_LengthGaugeFem, MCTDHX with one orbital (1s start) and
!> method.mb.orbBased.exactFieldFreeQ, propagator.splitStep.order2 (strang)
!> whose linear integrator is SIL in imaginary time (sil1.imagTimeQ): the
!> first step starts with the exact part, called with imagTimeQ_, and must
!> stop there instead of applying a real-time Cayley step.
!> ctest-pass-regex: the Cayley exact part .non-Ylm grids. is real-time only
program T_Method_FieldFreeCayleyImagTimeFem
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator

  implicit none

  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)
  real(R64) :: dt

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Method/T_Method_FieldFreeCayleyImagTimeFem.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
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
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup

  call Say_Fabricate("program")
  dt = Json_Get("program.timeStep", 0.01_R64)

  !==========================================
  call Say_Section("perform action")
  !==========================================

  call Propagator_Propagate(Method_state, 0.0_R64, dt)

  error stop "T_Method_FieldFreeCayleyImagTimeFem failure: imaginary time was accepted by the Cayley exact part"

end program
