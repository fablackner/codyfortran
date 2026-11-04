program T_Grid_Fem_Sliver
  !-----------------------------------------------------------------------------
  ! Near-duplicate centers: C1 = (0.3, -0.2, 0.5) (hMin 0.25) and
  ! C2 = (0.3000001, -0.2000001, 0.5000001) (hMin 0.2) are 1e-7 apart along
  ! every axis. The mesh grading merges them into C2 (the smaller hMin), so
  ! the mesh equals the one-center mesh of C2 (nLocals 5, box [-6, 6] x
  ! [-6, 6] x [-6, 6.5], growth 1.6, hMax 2 -> n = 47 per axis, N = 103823)
  ! and the Poisson solves stay accurate. Without the merge the 1e-7 sliver
  ! elements made the FDM spectrum negative (smallest eigenvalue -0.18 along
  ! z) and the Poisson solution garbage.
  !  (a) the element boundaries of all axes equal those of C2 alone, no
  !      element below 0.1 hMin, positive FDM spectra
  !  (b) default multipole center: the centroid of C1 and C2
  !  (c) Gaussian charge at C1: finite potential against erf(r)/r, self-energy
  !  (d) point charges -1 at C1 (1e-7 off the node: Lagrange basis values)
  !      and C2: finite potential, far field against -1/|r - C1| - 1/|r - C2|,
  !      value at the C2 node
  ! Prototype references (numpy, one-center mesh): (c) 2.1033e-5, self-energy
  ! error -4.566e-7; (d) 6.3737e-5 (r > 1.5), 1.1015e-5 (r > 3),
  ! -317.4305561012 at the C2 node (-317.4320509163 for charge -2 at C2).
  !-----------------------------------------------------------------------------
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: center1(3) = [0.3_R64, -0.2_R64, 0.5_R64]
  real(R64), parameter :: center2(3) = [0.3000001_R64, -0.2000001_R64, 0.5000001_R64]
  real(R64), parameter :: hMin2 = 0.2_R64

  integer(I32) :: iAxis, iNode
  real(R64) :: diffBounds, smallestElement, smallestEigenvalue, errGauss, selfEnergy, errPoint15, errPoint3
  real(R64) :: potNucleus
  real(R64), allocatable :: dist1(:), dist2(:)
  type(T_Grid_Fem_Axis) :: oneCenterAxes(3)
  complex(R64), allocatable :: load(:), pot(:)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  ! Reference: the mesh of C2 alone (element boundaries only)
  call Json_LoadJsonFile(manualFileName_="test/componentTests/Grid/Fem/T_Grid_Fem_Sliver_OneCenter.json", &
                         relativeToProjectDirQ_=.true.)
  call Grid_Fabricate
  oneCenterAxes = Grid_Fem_axes

  ! The near-duplicate pair
  call Json_LoadJsonFile(manualFileName_="test/componentTests/Grid/Fem/T_Grid_Fem_Sliver.json", &
                         relativeToProjectDirQ_=.true.)
  call Grid_Fabricate
  call Grid_Setup

  allocate (pot(Grid_nPoints))

  !==========================================
  call Say_Section("perform action")
  !==========================================

  !------------------------------------
  ! (a) the merged mesh is the one-center mesh of C2: identical element
  ! boundaries, no sliver element, positive FDM spectra
  !------------------------------------

  diffBounds = 0.0_R64
  smallestElement = huge(1.0_R64)
  smallestEigenvalue = huge(1.0_R64)
  do iAxis = 1, 3
    associate (axis => Grid_Fem_axes(iAxis))
      call check(error, size(axis % bounds), size(oneCenterAxes(iAxis) % bounds))
      if (allocated(error)) error stop "T_Grid_Fem_Sliver element count failure"
      diffBounds = max(diffBounds, maxval(abs(axis % bounds - oneCenterAxes(iAxis) % bounds)))
      smallestElement = min(smallestElement, minval(axis % bounds(1:) - axis % bounds(:axis % nElements - 1)))
      smallestEigenvalue = min(smallestEigenvalue, axis % eigenvalues(1))
    end associate
  end do

  call check(error, diffBounds < tiny(1.0_R64), "merged mesh differs from the one-center mesh")
  if (allocated(error)) error stop "T_Grid_Fem_Sliver merged mesh failure"
  call check(error, smallestElement > 0.1_R64 * hMin2, "sliver element")
  if (allocated(error)) error stop "T_Grid_Fem_Sliver smallest element failure"
  call check(error, smallestEigenvalue > 0.0_R64, "non-positive FDM eigenvalue")
  if (allocated(error)) error stop "T_Grid_Fem_Sliver FDM spectrum failure"

  !------------------------------------
  ! (b) default expansion center of the multipole boundary values: the
  ! centroid of the centers
  !------------------------------------

  call check(error, maxval(abs(Grid_Fem_multipoleCenter - 0.5_R64 * (center1 + center2))), 0.0_R64, thr=1.0e-15_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Sliver default multipole center failure"

  !------------------------------------
  ! (c) unit Gaussian charge at C1: V = erf(r)/r, self-energy 1/sqrt(2 pi)
  !------------------------------------

  dist1 = Distances(center1)
  load = Grid_Fem_weights * PI**(-1.5_R64) * exp(-dist1**2)
  call Grid_Fem_SolvePoisson(pot, load)

  call check(error, all(ieee_is_finite(real(pot, kind=R64))), "non-finite potential")
  if (allocated(error)) error stop "T_Grid_Fem_Sliver Gaussian potential finite failure"

  errGauss = maxval(abs(real(pot, kind=R64) - ErfOverR(dist1)))
  selfEnergy = 0.5_R64 * sum(real(load, kind=R64) * real(pot, kind=R64))

  call check(error, errGauss, 0.0_R64, thr=6.5e-5_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Sliver Gaussian potential failure"
  call check(error, selfEnergy, 1.0_R64 / sqrt(2.0_R64 * PI), thr=1.4e-6_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Sliver self-energy failure"

  !------------------------------------
  ! (d) point charges -1 at C1 (off the node by 1e-7) and C2 (node)
  !------------------------------------

  call Grid_Fem_SolvePoissonPointCharges(pot, [-1.0_R64, -1.0_R64], reshape([center1, center2], [3, 2]))

  call check(error, all(ieee_is_finite(real(pot, kind=R64))), "non-finite potential")
  if (allocated(error)) error stop "T_Grid_Fem_Sliver point charges finite failure"

  dist2 = Distances(center2)
  errPoint15 = maxval(abs(real(pot, kind=R64) + 1.0_R64 / max(dist1, 1.0e-12_R64) + 1.0_R64 / max(dist2, 1.0e-12_R64)), &
                      mask=dist2 > 1.5_R64)
  errPoint3 = maxval(abs(real(pot, kind=R64) + 1.0_R64 / max(dist1, 1.0e-12_R64) + 1.0_R64 / max(dist2, 1.0e-12_R64)), &
                     mask=dist2 > 3.0_R64)
  iNode = minloc(dist2, 1)
  potNucleus = real(pot(iNode), kind=R64)

  call check(error, dist2(iNode) < tiny(1.0_R64), "C2 is not a node")
  if (allocated(error)) error stop "T_Grid_Fem_Sliver C2 node failure"
  call check(error, errPoint15, 0.0_R64, thr=2.0e-4_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Sliver point charges (r > 1.5) failure"
  call check(error, errPoint3, 0.0_R64, thr=3.5e-5_R64)
  if (allocated(error)) error stop "T_Grid_Fem_Sliver point charges (r > 3) failure"
  call check(error, potNucleus, -317.4305561012_R64, thr=1.0e-10_R64, rel=.true.)
  if (allocated(error)) error stop "T_Grid_Fem_Sliver point charges value at the C2 node failure"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(A, 3I4, A, I0)') " (a) n = ", Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, ", N = ", Grid_nPoints
  write (*, '(A, ES12.4)') " (a) max |bounds - one-center bounds|          = ", diffBounds
  write (*, '(A, F12.6)') " (a) smallest element                          = ", smallestElement
  write (*, '(A, F12.6)') " (a) smallest FDM eigenvalue                   = ", smallestEigenvalue
  write (*, '(A, 3F12.8)') " (b) multipole center                          = ", Grid_Fem_multipoleCenter
  write (*, '(A, ES12.4)') " (c) Gaussian charge, max |V - erf(r)/r|       = ", errGauss
  write (*, '(A, F16.12, A, ES11.3, A)') " (c) self-energy = ", selfEnergy, "  (error ", &
    selfEnergy - 1.0_R64 / sqrt(2.0_R64 * PI), ")"
  write (*, '(A, ES12.4)') " (d) point charges, max |V - V_exact| (r > 1.5) = ", errPoint15
  write (*, '(A, ES12.4)') " (d) point charges, max |V - V_exact| (r > 3)   = ", errPoint3
  write (*, '(A, F16.10)') " (d) potential at the C2 node                   = ", potNucleus
  print *
  print *, "all fem sliver checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Distances |r_i - r0| of all nodes.
  function Distances(r0) result(res)
    real(R64), intent(in) :: r0(3)
    real(R64), allocatable :: res(:)

    integer(I32) :: jx, jy, jz

    allocate (res(Grid_nPoints))
    do jz = 1, Grid_Fem_nz
      do jy = 1, Grid_Fem_ny
        do jx = 1, Grid_Fem_nx
          res(Grid_Fem_Index(jx, jy, jz)) = norm2([Grid_Fem_axes(1) % points(jx), Grid_Fem_axes(2) % points(jy), &
                                                   Grid_Fem_axes(3) % points(jz)] - r0)
        end do
      end do
    end do

  end function

  !> Potential erf(r)/r of the unit Gaussian charge (2/sqrt(pi) at r = 0).
  elemental function ErfOverR(r) result(v)
    real(R64), intent(in) :: r
    real(R64) :: v

    if (r > 1.0e-12_R64) then
      v = erf(r) / r
    else
      v = 2.0_R64 / sqrt(PI)
    end if

  end function

end program
