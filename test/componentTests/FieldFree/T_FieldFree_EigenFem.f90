!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Expected error: fieldFree.eigen on grid.fem. The eigen variant diagonalizes
!> the radial blocks of h0 per l, which only Ylm grids have; the spectral cut
!> (fieldFree.eigen.cut) projects onto these eigenvectors, so it exists for
!> this variant alone. FieldFree_Eigen_Fabricate must stop on any other grid
!> (H atom on the mesh of T_Method_Sb_LengthGaugeFem, fieldFree.eigen with
!> cut 5).
!> ctest-pass-regex: fieldFree.eigen requires grid.ylm
program T_FieldFree_EigenFem
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_FieldFree

  implicit none

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/FieldFree/T_FieldFree_EigenFem.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call FieldFree_Fabricate

  error stop "T_FieldFree_EigenFem failure: fieldFree.eigen was accepted on grid.fem"

end program
