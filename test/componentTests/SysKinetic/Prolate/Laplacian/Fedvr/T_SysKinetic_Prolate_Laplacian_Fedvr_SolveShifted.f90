!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Exact shifted kinetic solve of the prolate FEDVR Laplacian
!> (SysKinetic_SolveShifted, (T + shift)⁻¹ by fast diagonalization of the
!> separated xi and eta operators): for a function with components in every
!> azimuthal channel (mmax = 2, graded xi elements), x = SolveShifted(f, shift)
!> must satisfy (T + shift) x = f to round-off, for a positive and a small
!> negative real shift and for a body type of another mass.
program T_SysKinetic_Prolate_Laplacian_Fedvr_SolveShifted
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Prolate
  use M_SysKinetic
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  complex(R64), allocatable :: f(:)
  real(R64) :: dev, xi, eta
  integer(I32) :: iGrid, m

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/SysKinetic/Prolate/Laplacian/Fedvr/T_SysKinetic_Prolate_Laplacian_Fedvr_SolveShifted.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate

  call Grid_Setup
  call SysKinetic_Setup

  call check(error, associated(SysKinetic_SolveShifted), "SysKinetic_SolveShifted must be bound")
  if (allocated(error)) error stop "T_SysKinetic_Prolate_Laplacian_Fedvr_SolveShifted failure: not bound"

  ! (xi^2 - 1)^(|m|/2) exp(-xi) (1 + eta/2 + i m/3) in every channel (zero on
  ! the axis for m /= 0, the channel Dirichlet point)
  allocate (f(Grid_nPoints))
  do iGrid = 1, Grid_nPoints
    xi = Grid_Prolate_xiCoord(iGrid)
    eta = Grid_Prolate_etaCoord(iGrid)
    m = Grid_Prolate_mCoord(iGrid)
    f(iGrid) = (xi**2 - 1.0_R64)**(0.5_R64 * abs(m)) * exp(-xi) * &
               cmplx(1.0_R64 + 0.5_R64 * eta, m / 3.0_R64, kind=R64)
  end do

  dev = Residual((0.7_R64, 0.0_R64), 1)
  print *, "real shift 0.7, body type 1: relative residual ", dev
  call check(error, dev < 1e-10_R64)
  if (allocated(error)) error stop "T_SysKinetic_Prolate_Laplacian_Fedvr_SolveShifted failure: real shift"

  dev = Residual((-0.05_R64, 0.0_R64), 1)
  print *, "real shift -0.05, body type 1: relative residual ", dev
  call check(error, dev < 1e-10_R64)
  if (allocated(error)) error stop "T_SysKinetic_Prolate_Laplacian_Fedvr_SolveShifted failure: negative shift"

  dev = Residual((0.7_R64, 0.0_R64), 2)
  print *, "real shift 0.7, body type 2 (mass 2): relative residual ", dev
  call check(error, dev < 1e-10_R64)
  if (allocated(error)) error stop "T_SysKinetic_Prolate_Laplacian_Fedvr_SolveShifted failure: body type 2"

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
