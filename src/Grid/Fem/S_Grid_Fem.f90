! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Fem) S_Grid_Fem
  !-----------------------------------------------------------------------------
  ! Mesh of the FEM grid: Grid_Fem_Fabricate (element partitions of the three
  ! axes, graded toward the centers or explicit), Grid_Setup (orchestration:
  ! Grid_Fem_SetupAxes, the 3D weights, the mesh summary, then
  ! Grid_Fem_SetupPoisson), the inner product and the node layout
  ! (Grid_Fem_Index, Grid_Fem_FillNodeCoordinate). The Laplacian and its fast
  ! diagonalization are in S_Grid_Fem_Laplacian, the Poisson solves in its
  ! child S_Grid_Fem_Poisson.
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Names of the x, y and z axes (keys grid.fem.<axis>Bounds, messages).
  character(len=1), parameter :: AXIS_NAMES(3) = ['x', 'y', 'z']

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fem_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid

    real(R64), allocatable :: centersX(:), centersY(:), centersZ(:)
    real(R64), allocatable :: bounds(:), multipoleCenter(:), defaultCenter(:)
    real(R64) :: growth, hMax, halfHMin
    integer(I32) :: iAxis, iCenter
    character(len=256) :: message

    call Say_Fabricate("grid.fem")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Grid_Fem_nLocals = Json_Get("grid.fem.nLocals", 6)
    if (Grid_Fem_nLocals < 2) error stop "grid.fem.nLocals must be >= 2"

    Grid_Fem_axes(1) % lowerBound = Json_Get("grid.fem.xmin", -10.0_R64)
    Grid_Fem_axes(1) % upperBound = Json_Get("grid.fem.xmax", 10.0_R64)
    Grid_Fem_axes(2) % lowerBound = Json_Get("grid.fem.ymin", -10.0_R64)
    Grid_Fem_axes(2) % upperBound = Json_Get("grid.fem.ymax", 10.0_R64)
    Grid_Fem_axes(3) % lowerBound = Json_Get("grid.fem.zmin", -10.0_R64)
    Grid_Fem_axes(3) % upperBound = Json_Get("grid.fem.zmax", 10.0_R64)

    do iAxis = 1, 3
      if (.not. (Grid_Fem_axes(iAxis) % upperBound > Grid_Fem_axes(iAxis) % lowerBound)) then
        error stop "grid.fem: the box needs xmin < xmax, ymin < ymax and zmin < zmax"
      end if
    end do

    ! Centers (nuclei) as parallel coordinate arrays, smallest element size at each
    centersX = Json_Get("grid.fem.centers.x", [real(R64) ::])
    centersY = Json_Get("grid.fem.centers.y", [real(R64) ::])
    centersZ = Json_Get("grid.fem.centers.z", [real(R64) ::])

    Grid_Fem_nCenters = size(centersX)
    if (size(centersY) .ne. Grid_Fem_nCenters .or. size(centersZ) .ne. Grid_Fem_nCenters) then
      error stop "grid.fem.centers: x, y and z must have equal sizes"
    end if

    if (allocated(Grid_Fem_centers)) deallocate (Grid_Fem_centers)
    allocate (Grid_Fem_centers(3, Grid_Fem_nCenters))
    Grid_Fem_centers(1, :) = centersX
    Grid_Fem_centers(2, :) = centersY
    Grid_Fem_centers(3, :) = centersZ

    Grid_Fem_centerHMin = Json_Get("grid.fem.centers.hMin", [(0.2_R64, iCenter=1, Grid_Fem_nCenters)])
    if (size(Grid_Fem_centerHMin) .ne. Grid_Fem_nCenters) error stop "grid.fem.centers.hMin needs one entry per center"
    if (any(.not. (Grid_Fem_centerHMin > 0.0_R64))) error stop "grid.fem.centers.hMin must be > 0"

    ! Every center at least 0.5 hMin inside the box faces: a center on or next
    ! to a face would put a sliver element (or the nucleus) at the Dirichlet
    ! surface (the comparison also rejects NaN coordinates)
    do iCenter = 1, Grid_Fem_nCenters
      halfHMin = 0.5_R64 * Grid_Fem_centerHMin(iCenter)
      do iAxis = 1, 3
        if (.not. (Grid_Fem_centers(iAxis, iCenter) - Grid_Fem_axes(iAxis) % lowerBound >= halfHMin .and. &
                   Grid_Fem_axes(iAxis) % upperBound - Grid_Fem_centers(iAxis, iCenter) >= halfHMin)) then
          write (message, '(A, I0, A)') "grid.fem.centers: center ", iCenter, " is closer than 0.5*hMin to a box face "// &
            "along "//AXIS_NAMES(iAxis)//" (every center must lie at least 0.5*hMin inside the box)"
          error stop trim(message)
        end if
      end do
    end do

    ! Mesh grading toward the centers
    growth = Json_Get("grid.fem.growth", 1.5_R64)
    if (.not. (growth >= 1.0_R64)) error stop "grid.fem.growth must be >= 1"

    hMax = Json_Get("grid.fem.hMax", 2.0_R64)
    if (.not. (hMax > 0.0_R64)) error stop "grid.fem.hMax must be > 0"

    ! Element partition of every axis: explicit bounds or graded toward the
    ! centers (generated here so that Grid_nPoints is known after fabrication)
    do iAxis = 1, 3
      associate (axis => Grid_Fem_axes(iAxis))

        if (Json_GetExistence("grid.fem."//AXIS_NAMES(iAxis)//"Bounds")) then
          bounds = Json_Get("grid.fem."//AXIS_NAMES(iAxis)//"Bounds", [real(R64) ::])
          call ValidateExplicitBounds(bounds, axis % lowerBound, axis % upperBound, "grid.fem."//AXIS_NAMES(iAxis)//"Bounds")

          ! Galerkin point charges work anywhere, but a nucleus off the element
          ! vertices loses the C0 resolution of its cusp
          do iCenter = 1, Grid_Fem_nCenters
            if (.not. any(abs(bounds - Grid_Fem_centers(iAxis, iCenter)) <= 0.0_R64)) then
              write (*, '(A, I0, A)') "  note: grid.fem."//AXIS_NAMES(iAxis)//"Bounds misses the coordinate of center ", &
                iCenter, " (not at an element vertex)"
            end if
          end do
        else
          call FillGradedBounds(bounds, axis % lowerBound, axis % upperBound, Grid_Fem_centers(iAxis, :), &
                                Grid_Fem_centerHMin, growth, hMax, AXIS_NAMES(iAxis))
        end if

        axis % nElements = size(bounds) - 1
        axis % nPoints = axis % nElements * (Grid_Fem_nLocals - 1) - 1
        if (axis % nPoints < 1) error stop "grid.fem: every axis needs an interior node (more elements or a larger nLocals)"

        if (allocated(axis % bounds)) deallocate (axis % bounds)
        allocate (axis % bounds(0:axis % nElements))
        axis % bounds(:) = bounds(:)

      end associate
    end do

    Grid_Fem_nx = Grid_Fem_axes(1) % nPoints
    Grid_Fem_ny = Grid_Fem_axes(2) % nPoints
    Grid_Fem_nz = Grid_Fem_axes(3) % nPoints
    Grid_nPoints = Grid_Fem_nx * Grid_Fem_ny * Grid_Fem_nz

    ! Multipole expansion of the free-space Poisson boundary values
    Grid_Fem_multipoleLmax = Json_Get("grid.fem.multipoleLmax", 6)
    if (Grid_Fem_multipoleLmax < 0) error stop "grid.fem.multipoleLmax must be >= 0"

    ! Expansion center: by default the centroid of the centers (nuclei), where
    ! the charge of a molecule sits; the box center if there are none
    if (Grid_Fem_nCenters > 0) then
      defaultCenter = sum(Grid_Fem_centers, dim=2) / real(Grid_Fem_nCenters, R64)
    else
      defaultCenter = 0.5_R64 * (Grid_Fem_axes(:) % lowerBound + Grid_Fem_axes(:) % upperBound)
    end if

    multipoleCenter = Json_Get("grid.fem.multipoleCenter", defaultCenter)
    if (size(multipoleCenter) .ne. 3) error stop "grid.fem.multipoleCenter needs three coordinates"
    do iAxis = 1, 3
      if (.not. (multipoleCenter(iAxis) > Grid_Fem_axes(iAxis) % lowerBound .and. &
                 multipoleCenter(iAxis) < Grid_Fem_axes(iAxis) % upperBound)) then
        error stop "grid.fem.multipoleCenter must lie strictly inside the box"
      end if
    end do
    Grid_Fem_multipoleCenter(:) = multipoleCenter(:)

    Grid_InnerProduct => InnerProduct
    Grid_Setup => Setup

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Nodes, weights and 1D operators of the three axes (Grid_Fem_SetupAxes),
  !> 3D weights, a mesh summary and the tables of the multipole boundary
  !> values (solid-harmonic polynomials, face-node geometry and irregular solid
  !> harmonics there, convergence radius; Grid_Fem_SetupPoisson, which ends the
  !> summary with the multipole lines).
  subroutine Setup
    use M_Utils_Say
    use M_Grid

    integer(I32) :: ix, iy, iz, iGrid

    call Say_Setup("grid.fem")

    call Grid_Fem_SetupAxes

    ! 3D product weights in the node layout (x fastest)
    if (allocated(Grid_Fem_weights)) deallocate (Grid_Fem_weights)
    allocate (Grid_Fem_weights(Grid_nPoints))

    iGrid = 0
    do iz = 1, Grid_Fem_nz
      do iy = 1, Grid_Fem_ny
        do ix = 1, Grid_Fem_nx
          iGrid = iGrid + 1
          Grid_Fem_weights(iGrid) = Grid_Fem_axes(1) % weights(ix) * Grid_Fem_axes(2) % weights(iy) * &
                                    Grid_Fem_axes(3) % weights(iz)
        end do
      end do
    end do

    call PrintMeshSummary
    call Grid_Fem_SetupPoisson

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Short mesh summary: elements, interior nodes and element sizes per axis
  !> (Grid_Fem_SetupPoisson adds the lines of the multipole boundary values).
  subroutine PrintMeshSummary
    use M_Grid

    real(R64), allocatable :: sizes(:)
    integer(I32) :: iAxis

    write (*, '(A, I0, A, I0)') "  grid.fem: nLocals = ", Grid_Fem_nLocals, ", multipoleLmax = ", Grid_Fem_multipoleLmax
    write (*, '(A)') "  axis       lower       upper  elements    nodes   smallest    largest"
    do iAxis = 1, 3
      associate (axis => Grid_Fem_axes(iAxis))
        sizes = axis % bounds(1:axis % nElements) - axis % bounds(0:axis % nElements - 1)
        write (*, '(5X, A1, 2F12.4, I10, I9, 2F11.4)') AXIS_NAMES(iAxis), axis % lowerBound, axis % upperBound, &
          axis % nElements, axis % nPoints, minval(sizes), maxval(sizes)
      end associate
    end do
    write (*, '(A, I0)') "  interior nodes in total: ", Grid_nPoints

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Explicit element boundaries of one axis (grid.fem.<axis>Bounds): at least
  !> two entries, strictly increasing, first and last on the box bounds (up to
  !> 1e-10 relative to the box length, then snapped onto them). The checks are
  !> NaN-safe (a NaN entry fails them).
  subroutine ValidateExplicitBounds(bounds, lower, upper, key)
    real(R64), intent(inout)     :: bounds(:)
    real(R64), intent(in)        :: lower
    real(R64), intent(in)        :: upper
    character(len=*), intent(in) :: key

    real(R64) :: tolerance
    integer(I32) :: n

    n = size(bounds)
    if (n < 2) error stop key//" needs at least two entries"

    tolerance = 1.0e-10_R64 * (upper - lower)
    if (.not. (abs(bounds(1) - lower) <= tolerance .and. abs(bounds(n) - upper) <= tolerance)) then
      error stop key//" must start and end at the box bounds"
    end if
    bounds(1) = lower
    bounds(n) = upper

    if (.not. all(bounds(2:n) > bounds(1:n - 1))) error stop key//" must be strictly increasing"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Element boundaries of one axis graded toward the center coordinates c_k,
  !> with the target element size
  !>   s(x) = min(hMax, min_k(h_k + (growth - 1)|x - c_k|)).
  !> Center coordinates closer than 0.25 min(h_k, h_l) are merged (each one
  !> against the kept coordinate of the current cluster, in ascending order):
  !> the cluster keeps the center with the smallest h_k (on ties the lowest
  !> index), and s(x) runs over the kept centers only, so a near duplicate
  !> gives exactly the mesh of one center instead of a sliver element; a note
  !> names every merged center off the kept coordinate (not at an element
  !> vertex; its Galerkin point charge stays exact through the basis values).
  !> Mandatory boundaries are the box bounds and the kept coordinates. A gap
  !> [p, q] between consecutive mandatory boundaries gets
  !> n = max(1, ceiling(F - 1e-9)) elements, F = int_p^q dx/s(x) (trapezoidal
  !> rule on 4000 subintervals), split where the cumulative integral reaches
  !> j F/n (linear inverse interpolation).
  subroutine FillGradedBounds(bounds, lower, upper, centers, hMins, growth, hMax, axisName)
    use stdlib_sorting, only: sort_index

    real(R64), allocatable, intent(out) :: bounds(:)
    real(R64), intent(in)               :: lower
    real(R64), intent(in)               :: upper
    real(R64), intent(in)               :: centers(:)
    real(R64), intent(in)               :: hMins(:)
    real(R64), intent(in)               :: growth
    real(R64), intent(in)               :: hMax
    character(len=*), intent(in)        :: axisName

    integer(I32), parameter :: nSub = 4000
    real(R64), parameter :: mergeFraction = 0.25_R64

    real(R64), allocatable :: sortedCenters(:), keptCenters(:), keptHMins(:), mandatory(:), xs(:), cumulative(:)
    real(R64) :: gapStart, gapEnd, step, fPrev, fCur, level, slope
    integer(I32), allocatable :: order(:), kept(:), clusterOf(:)
    integer(I32) :: nCenters, nKept, k, kNew, kKept, iGap, i, j, nElements

    nCenters = size(centers)

    ! Centers in ascending coordinate order (stable: ties keep the index order)
    sortedCenters = centers
    allocate (order(nCenters), kept(nCenters), clusterOf(nCenters))
    if (nCenters > 0) call sort_index(sortedCenters, order)

    ! Clusters of near duplicates and their kept centers
    nKept = 0
    do k = 1, nCenters
      kNew = order(k)

      if (nKept > 0) then
        kKept = kept(nKept)
        if (centers(kNew) - centers(kKept) < mergeFraction * min(hMins(kNew), hMins(kKept))) then
          clusterOf(kNew) = nKept
          ! keep the smaller hMin, on ties the lower index
          if (hMins(kNew) < hMins(kKept) .or. (.not. hMins(kNew) > hMins(kKept) .and. kNew < kKept)) kept(nKept) = kNew
          cycle
        end if
      end if

      nKept = nKept + 1
      kept(nKept) = kNew
      clusterOf(kNew) = nKept
    end do

    do k = 1, nCenters
      kKept = kept(clusterOf(k))
      if (abs(centers(k) - centers(kKept)) > 0.0_R64) then
        write (*, '(A, I0, A, I0, A, ES9.2, A)') "  note: grid.fem.centers: along "//axisName//" center ", k, &
          " is merged with center ", kKept, " (distance ", abs(centers(k) - centers(kKept)), &
          "), so it is not at an element vertex (its Galerkin point charge stays exact through the basis values)"
      end if
    end do

    keptCenters = centers(kept(1:nKept))
    keptHMins = hMins(kept(1:nKept))

    ! Mandatory boundaries: box bounds and kept center coordinates
    mandatory = [lower, keptCenters, upper]

    allocate (xs(0:nSub), cumulative(0:nSub))
    bounds = [lower]

    do iGap = 1, size(mandatory) - 1
      gapStart = mandatory(iGap)
      gapEnd = mandatory(iGap + 1)

      ! Cumulative integral of 1/s on nSub equal subintervals
      step = (gapEnd - gapStart) / nSub
      do i = 0, nSub - 1
        xs(i) = gapStart + i * step
      end do
      xs(nSub) = gapEnd

      cumulative(0) = 0.0_R64
      fPrev = 1.0_R64 / TargetSize(xs(0))
      do i = 1, nSub
        fCur = 1.0_R64 / TargetSize(xs(i))
        cumulative(i) = cumulative(i - 1) + 0.5_R64 * (fCur + fPrev) * (xs(i) - xs(i - 1))
        fPrev = fCur
      end do

      nElements = max(1, ceiling(cumulative(nSub) - 1.0e-9_R64))

      ! Interior boundaries where the cumulative integral reaches j F/n
      i = 1
      do j = 1, nElements - 1
        level = j * (cumulative(nSub) / nElements)
        do while (cumulative(i) < level)
          i = i + 1
        end do
        slope = (xs(i) - xs(i - 1)) / (cumulative(i) - cumulative(i - 1))
        bounds = [bounds, slope * (level - cumulative(i - 1)) + xs(i - 1)]
      end do

      bounds = [bounds, gapEnd]
    end do

  contains

    !> Target element size s(x) (kept centers)
    pure function TargetSize(x) result(s)
      real(R64), intent(in) :: x
      real(R64) :: s

      integer(I32) :: kCenter

      s = hMax
      do kCenter = 1, size(keptCenters)
        s = min(s, keptHMins(kCenter) + (growth - 1.0_R64) * abs(x - keptCenters(kCenter)))
      end do

    end function

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Inner product with the 3D product weights (plain accumulation, no
  !> temporaries).
  pure function InnerProduct(fConjg, f) result(res)

    complex(R64), intent(in), contiguous :: fConjg(:)
    complex(R64), intent(in), contiguous :: f(:)
    complex(R64) :: res

    integer(I32) :: i

    res = (0.0_R64, 0.0_R64)
    do i = 1, size(f)
      res = res + conjg(fConjg(i)) * f(i) * Grid_Fem_weights(i)
    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure module function Grid_Fem_Index(ix, iy, iz) result(i)

    integer(I32), intent(in) :: ix
    integer(I32), intent(in) :: iy
    integer(I32), intent(in) :: iz
    integer(I32) :: i

    i = ix + Grid_Fem_nx * (iy - 1) + Grid_Fem_nx * Grid_Fem_ny * (iz - 1)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fem_FillNodeCoordinate(coord, direction)
    !$ use omp_lib, only: omp_in_parallel
    use M_Grid

    real(R64), allocatable, intent(out) :: coord(:)
    real(R64), intent(in)               :: direction(3)

    real(R64) :: yz
    integer(I32) :: ix, iy, iz

    if (.not. allocated(Grid_Fem_axes(3) % points)) then
      error stop "Grid_Fem_FillNodeCoordinate: call Grid_Setup first"
    end if

    allocate (coord(Grid_nPoints))

    !$omp parallel do default(shared) private(ix, iy, iz, yz) schedule(static) if (.not. omp_in_parallel())
    do iz = 1, Grid_Fem_nz
      do iy = 1, Grid_Fem_ny
        yz = direction(2) * Grid_Fem_axes(2) % points(iy) + direction(3) * Grid_Fem_axes(3) % points(iz)
        do ix = 1, Grid_Fem_nx
          coord(Grid_Fem_Index(ix, iy, iz)) = direction(1) * Grid_Fem_axes(1) % points(ix) + yz
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine

end submodule
