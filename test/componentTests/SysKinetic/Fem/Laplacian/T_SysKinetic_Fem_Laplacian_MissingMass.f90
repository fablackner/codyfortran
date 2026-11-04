program T_SysKinetic_Fem_Laplacian_MissingMass
  !-----------------------------------------------------------------------------
  ! sysKinetic.fem.laplacian with several masses ([1.0, 2.5]) must cover every
  ! body type a caller requests: T for body type 3 stops with a clear message
  ! (a single bodyMass entry would serve every body type, see
  ! T_SysKinetic_Fem_Laplacian).
  !-----------------------------------------------------------------------------
!> ctest-pass-regex: no mass for body type 3
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic

  implicit none

  complex(R64), allocatable :: f(:), tf(:)

  call Say_Hello

  call Json_LoadJsonFile(manualFileName_= &
                         "test/componentTests/SysKinetic/Fem/Laplacian/T_SysKinetic_Fem_Laplacian_MissingMass.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate

  call Grid_Setup
  call SysKinetic_Setup

  allocate (f(Grid_nPoints), tf(Grid_nPoints))
  f(:) = (1.0_R64, 0.0_R64)

  ! must stop with the expected message
  call SysKinetic_MultiplyWithKineticOp(tf, f, 0.0_R64, bt_=3)

  print *, "T_SysKinetic_Fem_Laplacian_MissingMass: body type 3 without a mass was not rejected"

end program
