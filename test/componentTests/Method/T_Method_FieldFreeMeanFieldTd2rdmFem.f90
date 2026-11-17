!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Expected error: fieldFree.cayley.meanFieldQ with TD-2RDM. The frozen
!> mean-field reference of the Cayley exact part is subtracted in the
!> nonlinear part of the methods that declare it
!> (Method_Mb_OrbBased_implementsMeanFieldReferenceQ: MCTDHX, TDKS); TD-2RDM
!> binds Method_PropagateExactPart as well, so without the check its exact part
!> would propagate h0 + U while its nonlinear part keeps the full interaction.
!> Method_Mb_OrbBased_Fabricate must stop (H atom on the mesh of
!> T_Method_Sb_LengthGaugeFem with a Coulomb interaction, td2rdm with
!> fieldFree.cayley.meanFieldQ).
!> ctest-pass-regex: fieldFree.cayley.meanFieldQ is not implemented for the selected method
program T_Method_FieldFreeMeanFieldTd2rdmFem
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_FieldFree
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
  call FieldFree_Fabricate
  call Method_Fabricate

  error stop "T_Method_FieldFreeMeanFieldTd2rdmFem failure: fieldFree.cayley.meanFieldQ was accepted with td2rdm"

end program
