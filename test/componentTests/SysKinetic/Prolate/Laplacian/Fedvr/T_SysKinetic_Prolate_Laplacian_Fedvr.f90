program T_SysKinetic_Prolate_Laplacian_Fedvr
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Prolate
  use M_SysKinetic
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: iGrid
  complex(R64) :: braKetF, braKetG, expectation
  complex(R64), allocatable :: fOrb(:), gOrb(:), tfOrb(:), tgOrb(:)
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/SysKinetic/Prolate/Laplacian/Fedvr/T_SysKinetic_Prolate_Laplacian_Fedvr.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate

  call Grid_Setup
  call SysKinetic_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  ! Deterministic pseudo-random test fields populating all m channels
  allocate (fOrb(Grid_nPoints), gOrb(Grid_nPoints), tfOrb(Grid_nPoints), tgOrb(Grid_nPoints))
  do iGrid = 1, Grid_nPoints
    fOrb(iGrid) = cmplx(sin(0.137_R64 * iGrid), cos(0.291_R64 * iGrid), kind=R64)
    gOrb(iGrid) = cmplx(cos(0.113_R64 * iGrid), sin(0.407_R64 * iGrid), kind=R64)
  end do

  call SysKinetic_MultiplyWithKineticOp(tfOrb, fOrb, 0.0_R64)
  call SysKinetic_MultiplyWithKineticOp(tgOrb, gOrb, 0.0_R64)

  !------------------------------------
  ! Hermiticity with respect to the metric-weighted inner product:
  ! <f|T g> = conjg(<g|T f>)
  !------------------------------------

  braKetF = Grid_InnerProduct(fOrb, tgOrb)
  braKetG = Grid_InnerProduct(gOrb, tfOrb)

  call check(error, abs(braKetF - conjg(braKetG)) / abs(braKetF), 0.0_R64, thr=1e-12_R64)
  if (allocated(error)) error stop "T_SysKinetic_Prolate Hermiticity failure"

  !------------------------------------
  ! positivity: <f|T f> real and positive
  !------------------------------------

  expectation = Grid_InnerProduct(fOrb, tfOrb)

  call check(error, aimag(expectation) / abs(expectation), 0.0_R64, thr=1e-12_R64)
  if (allocated(error)) error stop "T_SysKinetic_Prolate real expectation failure"

  call check(error, real(expectation, kind=R64) > 0.0_R64)
  if (allocated(error)) error stop "T_SysKinetic_Prolate positivity failure"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "all prolate kinetic checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
