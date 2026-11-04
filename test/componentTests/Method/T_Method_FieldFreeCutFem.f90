!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Expected error: method.mb.orbBased.exactFieldFreeCut on grid.fem. The
!> spectral cut projects onto h0 eigenvectors, which only the eigen-decomposed
!> exact part of Ylm grids has; the Cayley exact part of grid.fem binds no
!> Method_Mb_OrbBased_ProjectBelowFieldFreeCut, so Method_Mb_OrbBased_Fabricate
!> must stop before it would bind Method_Mb_OrbBased_ApplySpectralCut to a
!> null pointer (H atom on the mesh of T_Method_Sb_LengthGaugeFem, MCTDHX with
!> exactFieldFreeQ and exactFieldFreeCut 5).
!> ctest-pass-regex: exactFieldFreeCut is only available on Ylm grids
program T_Method_FieldFreeCutFem
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_Method

  implicit none

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Method/T_Method_FieldFreeCutFem.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call Method_Fabricate

  error stop "T_Method_FieldFreeCutFem failure: exactFieldFreeCut was accepted on grid.fem"

end program
