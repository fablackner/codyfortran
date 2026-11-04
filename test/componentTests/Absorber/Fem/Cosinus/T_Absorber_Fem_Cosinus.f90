!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Component test for the boundary-layer absorber absorber.fem.cosinus on an
!> asymmetric box [-6, 6] x [-5, 7] x [-4, 4] (nLocals 5, hMax 1), width 2,
!> order 4.
!>
!> Applies the absorber to constant orbitals, which then equal the mask, and
!> checks:
!> - mask = 1 (exactly) at every node at least `width` away from all faces
!> - mask = prod_a cos^(1/order)(pi/2 (width - d_a)/width) (d_a < width, d_a
!>   the distance to the nearer face of axis a) at every node, in particular
!>   at boundary-layer nodes and in the edges/corners of the box
!> - monotone decay toward each face along the axis lines through the inner
!>   box, and 0 < mask <= 1
program T_Absorber_Fem_Cosinus
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use M_Absorber
  use M_Absorber_Fem_Cosinus
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: width = 2.0_R64
  integer(I32), parameter :: order = 4

  integer(I32) :: ix, iy, iz, i, iAxis, nInner, nLayer, iMid(3)
  real(R64) :: maskExpected, errMask, errInner, mask, profileRef, layerNodeMask, layerNodeX
  real(R64), allocatable :: line(:)
  logical :: monotoneQ, rangeQ
  complex(R64), allocatable :: orbs(:, :)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Absorber/Fem/Cosinus/T_Absorber_Fem_Cosinus.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Absorber_Fabricate

  call Grid_Setup
  call Absorber_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  call check(error, Absorber_Fem_Cosinus_width, width, thr=1.0e-15_R64)
  if (allocated(error)) error stop "T_Absorber_Fem_Cosinus width failure"
  call check(error, Absorber_Fem_Cosinus_order, order)
  if (allocated(error)) error stop "T_Absorber_Fem_Cosinus order failure"

  ! constant orbitals: afterwards they equal the mask (times the constant)
  allocate (orbs(Grid_nPoints, 2))
  orbs(:, 1) = (1.0_R64, 0.0_R64)
  orbs(:, 2) = (0.0_R64, 2.0_R64)

  call Absorber_ApplyAbsorber(orbs)

  !------------------------------------
  ! mask against the product formula; exactly 1 in the inner box
  !------------------------------------

  errMask = 0.0_R64
  errInner = 0.0_R64
  nInner = 0
  nLayer = 0
  rangeQ = .true.
  layerNodeMask = 1.0_R64
  layerNodeX = 0.0_R64

  do iz = 1, Grid_Fem_nz
    do iy = 1, Grid_Fem_ny
      do ix = 1, Grid_Fem_nx
        i = Grid_Fem_Index(ix, iy, iz)
        maskExpected = Profile(1, ix) * Profile(2, iy) * Profile(3, iz)
        mask = real(orbs(i, 1), kind=R64)

        errMask = max(errMask, abs(mask - maskExpected), abs(aimag(orbs(i, 2)) - 2.0_R64 * maskExpected))
        rangeQ = rangeQ .and. mask > 0.0_R64 .and. mask <= 1.0_R64

        if (maskExpected < 1.0_R64) then
          nLayer = nLayer + 1
        else
          nInner = nInner + 1
          errInner = max(errInner, abs(mask - 1.0_R64))
        end if
      end do
    end do
  end do

  call check(error, errMask, 0.0_R64, thr=1.0e-15_R64)
  if (allocated(error)) error stop "T_Absorber_Fem_Cosinus mask formula failure"
  call check(error, nInner > 0 .and. nLayer > 0 .and. errInner <= 0.0_R64, "the inner box must have mask exactly 1")
  if (allocated(error)) error stop "T_Absorber_Fem_Cosinus inner box failure"
  call check(error, rangeQ, "the mask must lie in (0, 1]")
  if (allocated(error)) error stop "T_Absorber_Fem_Cosinus mask range failure"

  ! one boundary-layer node explicitly: the last x node on the axis line
  ! through the center (0, 0.5, 0)
  iMid = [minloc(abs(Grid_Fem_axes(1) % points), dim=1), minloc(abs(Grid_Fem_axes(2) % points - 0.5_R64), dim=1), &
          minloc(abs(Grid_Fem_axes(3) % points), dim=1)]
  layerNodeX = Grid_Fem_axes(1) % points(Grid_Fem_nx)
  layerNodeMask = real(orbs(Grid_Fem_Index(Grid_Fem_nx, iMid(2), iMid(3)), 1), kind=R64)
  profileRef = cos(0.5_R64 * PI * (width - (6.0_R64 - layerNodeX)) / width)**(1.0_R64 / real(order, R64))

  call check(error, layerNodeMask, profileRef, thr=1.0e-15_R64)
  if (allocated(error)) error stop "T_Absorber_Fem_Cosinus boundary-layer node failure"

  !------------------------------------
  ! monotone decay toward the faces along the axis lines through iMid
  !------------------------------------

  monotoneQ = .true.
  do iAxis = 1, 3
    line = AxisLine(iAxis)
    ! non-decreasing up to the middle node, non-increasing after it
    do i = 2, size(line)
      if (i <= iMid(iAxis)) then
        monotoneQ = monotoneQ .and. line(i) >= line(i - 1)
      else
        monotoneQ = monotoneQ .and. line(i) <= line(i - 1)
      end if
    end do
    monotoneQ = monotoneQ .and. line(1) < 1.0_R64 .and. line(size(line)) < 1.0_R64
  end do

  call check(error, monotoneQ, "the mask must decay monotonically toward every face")
  if (allocated(error)) error stop "T_Absorber_Fem_Cosinus monotone decay failure"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(A, I0, A, I0, A, I0)') " nodes per axis                          = ", Grid_Fem_nx, ", ", Grid_Fem_ny, ", ", &
    Grid_Fem_nz
  write (*, '(A, I0, A, I0)') " inner / boundary-layer nodes            = ", nInner, " / ", nLayer
  write (*, '(A, ES12.4)') " max |mask - product formula|            = ", errMask
  write (*, '(A, F12.8, A, F12.8)') " mask at the last x node (x = ", layerNodeX, ")  = ", layerNodeMask
  print *
  print *, "all fem absorber checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Layer profile of axis iAxis at its node j.
  function Profile(iAxis, j) result(res)
    integer(I32), intent(in) :: iAxis
    integer(I32), intent(in) :: j
    real(R64) :: res

    real(R64) :: d

    d = min(Grid_Fem_axes(iAxis) % points(j) - Grid_Fem_axes(iAxis) % lowerBound, &
            Grid_Fem_axes(iAxis) % upperBound - Grid_Fem_axes(iAxis) % points(j))

    if (d < width) then
      res = cos(0.5_R64 * PI * (width - d) / width)**(1.0_R64 / real(order, R64))
    else
      res = 1.0_R64
    end if

  end function

  !> Mask values along axis iAxis through the node iMid.
  function AxisLine(iAxis) result(values)
    integer(I32), intent(in) :: iAxis
    real(R64), allocatable :: values(:)

    integer(I32) :: j, idx(3)

    allocate (values(Grid_Fem_axes(iAxis) % nPoints))
    do j = 1, Grid_Fem_axes(iAxis) % nPoints
      idx = iMid
      idx(iAxis) = j
      values(j) = real(orbs(Grid_Fem_Index(idx(1), idx(2), idx(3)), 1), kind=R64)
    end do

  end function

end program
