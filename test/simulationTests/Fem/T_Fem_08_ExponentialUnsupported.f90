!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Expected error: propagator.exponential with method.mb.orbBased.exactFieldFreeQ
!> on grid.fem. The exact part of non-Ylm grids is the Cayley step, which has
!> no phi functions, so MCTDHX binds no Method_ApplyExactPartPhi and
!> Propagator_Exponential_Fabricate stops with a message pointing to
!> propagator.splitStep (the configuration of T_Fem_05_H2FieldFree with the
!> exponential Krogstad propagator, fabricated in the order of the other
!> real-time tests; it stops before any setup).
!> ctest-pass-regex: exactFieldFreeQ on non-Ylm grids .e.g. grid.fem. is a Cayley step without phi functions
program T_Fem_08_ExponentialUnsupported
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

  implicit none

  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/simulationTests/Fem/T_Fem_08_ExponentialUnsupported.json", &
                         relativeToProjectDirQ_=.true.)

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
  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  error stop "T_Fem_08_ExponentialUnsupported failure: propagator.exponential was accepted with the Cayley exact part"

end program
