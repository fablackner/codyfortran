!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Exact shifted kinetic solve of the FEDVR Ylm Laplacian
!> (SysKinetic_SolveShifted, (T + shift)⁻¹ by per-l eigen-decomposition):
!> for a function with components in every (l, m) channel (lmax = mmax = 2),
!> x = SolveShifted(f, shift) must satisfy (T + shift) x = f to round-off,
!> for a real (preconditioner) and a complex (Crank–Nicolson-like) shift and
!> for a body type of another mass.
program T_SysKinetic_Ylm_Laplacian_Fedvr_SolveShifted
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Ylm
  use M_SysKinetic
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  complex(R64), allocatable :: f(:)
  real(R64) :: dev
  integer(I32) :: iGrid

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/SysKinetic/Ylm/Laplacian/Fedvr/T_SysKinetic_Ylm_Laplacian_Fedvr_SolveShifted.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate

  call Grid_Setup
  call SysKinetic_Setup

  call check(error, associated(SysKinetic_SolveShifted), "SysKinetic_SolveShifted must be bound")
  if (allocated(error)) error stop "T_SysKinetic_Ylm_Laplacian_Fedvr_SolveShifted failure: not bound"

  ! r^l exp(-r) (1 + i m/3 + l/5) in every channel
  allocate (f(Grid_nPoints))
  do iGrid = 1, Grid_nPoints
    f(iGrid) = Grid_Ylm_rCoord(iGrid)**Grid_Ylm_lCoord(iGrid) * exp(-Grid_Ylm_rCoord(iGrid)) * &
               cmplx(1.0_R64 + Grid_Ylm_lCoord(iGrid) / 5.0_R64, Grid_Ylm_mCoord(iGrid) / 3.0_R64, kind=R64)
  end do

  dev = Residual((0.7_R64, 0.0_R64), 1)
  print *, "real shift 0.7, body type 1: relative residual ", dev
  call check(error, dev < 1e-10_R64)
  if (allocated(error)) error stop "T_SysKinetic_Ylm_Laplacian_Fedvr_SolveShifted failure: real shift"

  dev = Residual((0.3_R64, -2.0_R64), 1)
  print *, "complex shift (0.3, -2), body type 1: relative residual ", dev
  call check(error, dev < 1e-10_R64)
  if (allocated(error)) error stop "T_SysKinetic_Ylm_Laplacian_Fedvr_SolveShifted failure: complex shift"

  dev = Residual((0.7_R64, 0.0_R64), 2)
  print *, "real shift 0.7, body type 2 (mass 2): relative residual ", dev
  call check(error, dev < 1e-10_R64)
  if (allocated(error)) error stop "T_SysKinetic_Ylm_Laplacian_Fedvr_SolveShifted failure: body type 2"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> max |(T + shift) x − f| / max |f| with x = SolveShifted(f, shift).
  function Residual(shift, bt) result(res)
    complex(R64), intent(in) :: shift
    integer(I32), intent(in) :: bt
    real(R64) :: res

    complex(R64), allocatable :: x(:), tx(:)

    allocate (x(Grid_nPoints), tx(Grid_nPoints))
    call SysKinetic_SolveShifted(x, f, shift, bt_=bt)
    call SysKinetic_MultiplyWithKineticOp(tx, x, 0.0_R64, bt_=bt)
    res = maxval(abs(tx + shift * x - f)) / maxval(abs(f))

  end function

end program
