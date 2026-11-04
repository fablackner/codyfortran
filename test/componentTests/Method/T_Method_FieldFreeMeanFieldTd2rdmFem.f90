!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Expected error: method.mb.orbBased.exactFieldFreeMeanFieldQ with TD-2RDM.
!> The frozen mean-field reference of the Cayley exact part is subtracted in
!> the nonlinear part of MCTDHX only (orbitals and CI coefficients); TD-2RDM
!> binds Method_PropagateExactPart as well, so without the check its exact part
!> would propagate h0 + U while its nonlinear part keeps the full interaction.
!> Method_Mb_OrbBased_Fabricate must stop (H atom on the mesh of
!> T_Method_Sb_LengthGaugeFem with a Coulomb interaction, td2rdm with
!> exactFieldFreeQ and exactFieldFreeMeanFieldQ).
!> ctest-pass-regex: exactFieldFreeMeanFieldQ is only implemented for method.mb.orbBased.mctdhx
program T_Method_FieldFreeMeanFieldTd2rdmFem
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Method

  implicit none

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Method/T_Method_FieldFreeMeanFieldTd2rdmFem.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate

  error stop "T_Method_FieldFreeMeanFieldTd2rdmFem failure: exactFieldFreeMeanFieldQ was accepted with td2rdm"

end program
