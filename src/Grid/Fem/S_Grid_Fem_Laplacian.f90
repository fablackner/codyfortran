! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Fem) S_Grid_Fem_Laplacian
  !-----------------------------------------------------------------------------
  ! Discrete Laplacian of the FEM grid and its fast diagonalization (FDM): the
  ! 1D operator tables of the axes (built by Grid_Fem_SetupAxes), the banded
  ! Laplacian, the eigenbasis transforms and the shifted solves, with the real
  ! kernels they share with the Poisson solves of the child submodule
  ! S_Grid_Fem_Poisson (mode products, spectral inverse, real and imaginary
  ! parts, argument checks). The child reaches the tables and kernels by host
  ! association; none of them is public.
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Names of the x, y and z axes (messages).
  character(len=1), parameter :: AXIS_NAMES(3) = ['x', 'y', 'z']

  !> 1D operator tables of one axis on its interior nodes (p = nLocals - 1).
  type :: T_AxisOperators
    !> FDM eigenvectors S (nPoints, nPoints): S^T K S = diag(eigenvalues),
    !> S^T W S = 1 (columns in the order of the eigenvalues of the axis).
    real(R64), allocatable :: eigenvectors(:, :)
    !> Banded W^-1 K: laplacianBand(d, i) = (W^-1 K)(i, i + d), d = -p..p,
    !> zero where i + d is outside 1..nPoints. The 1D nodal second derivative
    !> is -sum_d laplacianBand(d, i) f(i + d).
    real(R64), allocatable :: laplacianBand(:, :)
    !> Stiffness coupling of the interior nodes to the lower/upper box-boundary
    !> node, K_full(i, boundary node) (nPoints; nonzero only in the first/last
    !> element).
    real(R64), allocatable :: lowerBoundaryColumn(:)
    real(R64), allocatable :: upperBoundaryColumn(:)
  end type

  !> Operator tables of the x, y and z axes (Grid_Fem_SetupAxes).
  type(T_AxisOperators) :: axisOperators(3)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fem_SetupAxes

    real(R64) :: hMinSmallest
    integer(I32) :: iAxis

    ! Reference size of the sliver-element warning: the smallest requested
    ! element size (no warning without centers)
    hMinSmallest = 0.0_R64
    if (Grid_Fem_nCenters > 0) hMinSmallest = minval(Grid_Fem_centerHMin)

    do iAxis = 1, 3
      call SetupAxis(axisOperators(iAxis), Grid_Fem_axes(iAxis), Grid_Fem_nLocals, AXIS_NAMES(iAxis), hMinSmallest)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Nodes, weights and eigenvalues (axis) and the 1D operator tables (ops) of
  !> one axis from its element boundaries:
  !>  - GLL nodes and assembled weights (FEDVR context on the bounds, both box
  !>    nodes included), of which the interior ones are kept;
  !>  - the assembled stiffness K_full(i, j) = int L_i' L_j' dx from the element
  !>    matrices Kref/h_e, Kref the reference stiffness on [-1/2, 1/2] (GLL
  !>    quadrature, exact for the degree 2p - 2 integrand);
  !>  - the band of W^-1 K_II and the boundary columns K_full(interior, box node);
  !>  - the FDM eigendecomposition W^-1/2 K_II W^-1/2 = Q diag(lambda) Q^T,
  !>    S = W^-1/2 Q, so S^T K_II S = diag(lambda) and S^T W S = 1; stops
  !>    unless the smallest eigenvalue is positive (-d^2/dx^2 with Dirichlet
  !>    conditions is positive definite: a broken eigendecomposition).
  !> Warns about sliver elements below 0.1 hMinSmallest (if > 0), which make
  !> the operators ill-conditioned.
  subroutine SetupAxis(ops, axis, nLocals, axisName, hMinSmallest)
    use M_Utils_Fedvr, only: T_Fedvr_Ctx, Fedvr_CreateCtxFromBounds
    use M_Utils_DerivativeFedvr, only: T_DerivativeFedvr_Ctx, DerivativeFedvr_CreateCtx
    use M_Utils_LapackLib, only: LapackLib_DiagonalizeGeneric

    type(T_AxisOperators), intent(out)   :: ops
    type(T_Grid_Fem_Axis), intent(inout) :: axis
    integer(I32), intent(in)             :: nLocals
    character(len=*), intent(in)         :: axisName
    real(R64), intent(in)                :: hMinSmallest

    type(T_Fedvr_Ctx) :: fedvrCtx
    type(T_DerivativeFedvr_Ctx) :: derivativeCtx
    real(R64), allocatable :: stiffnessRef(:, :), stiffness(:, :), a(:, :), q(:, :)
    real(R64) :: smallestElement
    integer(I32) :: p, n, nFull, iE, i0, i, j, d

    p = nLocals - 1

    smallestElement = minval(axis % bounds(1:axis % nElements) - axis % bounds(0:axis % nElements - 1))
    if (smallestElement < 0.1_R64 * hMinSmallest) then
      write (*, '(A, ES10.3, A, ES10.3, A)') "  warning: grid.fem: the smallest "//axisName//" element (", smallestElement, &
        ") is below 0.1 min(centers.hMin) = ", 0.1_R64 * hMinSmallest, " (sliver element: ill-conditioned operators)"
    end if

    call Fedvr_CreateCtxFromBounds(fedvrCtx, axis % bounds, nLocals, .false., .false.)
    call DerivativeFedvr_CreateCtx(derivativeCtx, fedvrCtx)

    nFull = fedvrCtx % nPoints
    n = nFull - 2

    ! Reference stiffness int L_i' L_j' on [-1/2, 1/2] (secondOrderMatrix is its
    ! negative), symmetrized against round-off
    stiffnessRef = -0.5_R64 * (derivativeCtx % secondOrderMatrix + transpose(derivativeCtx % secondOrderMatrix))

    ! Assembled stiffness on all nodes (box nodes 1 and nFull included)
    allocate (stiffness(nFull, nFull))
    stiffness = 0.0_R64
    do iE = 1, axis % nElements
      i0 = fedvrCtx % elements(iE) % iStart - 1
      stiffness(i0 + 1:i0 + nLocals, i0 + 1:i0 + nLocals) = stiffness(i0 + 1:i0 + nLocals, i0 + 1:i0 + nLocals) + &
                                                            stiffnessRef / fedvrCtx % elements(iE) % size
    end do

    ! Interior nodes; the two box nodes carry the Dirichlet values
    axis % points = fedvrCtx % points(2:nFull - 1)
    axis % weights = fedvrCtx % weights(2:nFull - 1)
    ops % lowerBoundaryColumn = stiffness(2:nFull - 1, 1)
    ops % upperBoundaryColumn = stiffness(2:nFull - 1, nFull)

    ! Banded W^-1 K_II
    allocate (ops % laplacianBand(-p:p, n))
    ops % laplacianBand = 0.0_R64
    do i = 1, n
      do d = max(-p, 1 - i), min(p, n - i)
        ops % laplacianBand(d, i) = stiffness(i + 1, i + d + 1) / axis % weights(i)
      end do
    end do

    ! FDM eigendecomposition (DSYEV, ascending eigenvalues)
    allocate (a(n, n))
    do j = 1, n
      do i = 1, n
        a(i, j) = stiffness(i + 1, j + 1) / sqrt(axis % weights(i) * axis % weights(j))
      end do
    end do

    call LapackLib_DiagonalizeGeneric(axis % eigenvalues, q, a, .true.)

    if (.not. (axis % eigenvalues(1) > 0.0_R64)) then
      error stop "grid.fem: the smallest FDM eigenvalue of the "//axisName//" axis is not positive "// &
        "(degenerate sliver elements? check grid.fem.centers or grid.fem."//axisName//"Bounds)"
    end if

    allocate (ops % eigenvectors(n, n))
    do j = 1, n
      ops % eigenvectors(:, j) = q(:, j) / sqrt(axis % weights(:))
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fem_ApplyLaplacian(lapF, f)

    complex(R64), intent(out), contiguous :: lapF(:)
    complex(R64), intent(in), contiguous  :: f(:)

    call CheckFieldSize(size(lapF), "Grid_Fem_ApplyLaplacian")
    call CheckFieldSize(size(f), "Grid_Fem_ApplyLaplacian")

    call ApplyBandedLaplacian(lapF, f, axisOperators(1) % laplacianBand, axisOperators(2) % laplacianBand, &
                              axisOperators(3) % laplacianBand, Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, Grid_Fem_nLocals - 1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> lapF = -sum_a (W_a^-1 K_a) f on the 3D node array: the x stencil runs along
  !> the contiguous index, the y and z stencils combine whole x rows. Only the
  !> band entries of nodes sharing an element are visited (BandRange).
  subroutine ApplyBandedLaplacian(lapF, f, bandX, bandY, bandZ, nx, ny, nz, p)
    !$ use omp_lib, only: omp_in_parallel

    integer(I32), intent(in)  :: nx
    integer(I32), intent(in)  :: ny
    integer(I32), intent(in)  :: nz
    integer(I32), intent(in)  :: p
    complex(R64), intent(out) :: lapF(nx, ny, nz)
    complex(R64), intent(in)  :: f(nx, ny, nz)
    real(R64), intent(in)     :: bandX(-p:p, nx)
    real(R64), intent(in)     :: bandY(-p:p, ny)
    real(R64), intent(in)     :: bandZ(-p:p, nz)

    integer(I32) :: ix, iy, iz, d, dLo, dHi
    complex(R64) :: acc

    !$omp parallel do default(shared) private(ix, iy, iz, d, dLo, dHi, acc) schedule(static) if (.not. omp_in_parallel())
    do iz = 1, nz
      do iy = 1, ny

        ! x stencil
        do ix = 1, nx
          call BandRange(dLo, dHi, ix, nx, p)
          acc = (0.0_R64, 0.0_R64)
          do d = dLo, dHi
            acc = acc + bandX(d, ix) * f(ix + d, iy, iz)
          end do
          lapF(ix, iy, iz) = -acc
        end do

        ! y stencil
        call BandRange(dLo, dHi, iy, ny, p)
        do d = dLo, dHi
          lapF(:, iy, iz) = lapF(:, iy, iz) - bandY(d, iy) * f(:, iy + d, iz)
        end do

        ! z stencil
        call BandRange(dLo, dHi, iz, nz, p)
        do d = dLo, dHi
          lapF(:, iy, iz) = lapF(:, iy, iz) - bandZ(d, iz) * f(:, iy, iz + d)
        end do

      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Offsets dLo..dHi of the band entries of interior node i (1..n) that can be
  !> nonzero: node i sits at local position k = mod(i, p) of its element and
  !> couples to the nodes i - k .. i - k + p of that element (k > 0), an
  !> element-boundary node (k = 0) to i - p .. i + p; clipped to 1..n.
  pure subroutine BandRange(dLo, dHi, i, n, p)

    integer(I32), intent(out) :: dLo
    integer(I32), intent(out) :: dHi
    integer(I32), intent(in)  :: i
    integer(I32), intent(in)  :: n
    integer(I32), intent(in)  :: p

    integer(I32) :: k

    k = mod(i, p)
    if (k .eq. 0) then
      dLo = -p
      dHi = p
    else
      dLo = -k
      dHi = p - k
    end if

    dLo = max(dLo, 1 - i)
    dHi = min(dHi, n - i)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fem_ToEigenbasis(c, f)
    use M_Grid

    complex(R64), intent(out), contiguous :: c(:)
    complex(R64), intent(in), contiguous  :: f(:)

    real(R64), allocatable :: parts(:, :), res(:)

    call CheckFieldSize(size(c), "Grid_Fem_ToEigenbasis")
    call CheckFieldSize(size(f), "Grid_Fem_ToEigenbasis")

    call ExtractParts(parts, f, weightedQ=.true.)

    allocate (res(Grid_nPoints))

    call ModeProducts(res, parts(:, 1), Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, transposeQ=.true.)
    call SetRealPart(c, res)

    if (size(parts, 2) > 1) then
      call ModeProducts(res, parts(:, 2), Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, transposeQ=.true.)
      call SetImagPart(c, res)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fem_FromEigenbasis(f, c)
    use M_Grid

    complex(R64), intent(out), contiguous :: f(:)
    complex(R64), intent(in), contiguous  :: c(:)

    real(R64), allocatable :: parts(:, :), res(:)

    call CheckFieldSize(size(f), "Grid_Fem_FromEigenbasis")
    call CheckFieldSize(size(c), "Grid_Fem_FromEigenbasis")

    call ExtractParts(parts, c, weightedQ=.false.)

    allocate (res(Grid_nPoints))

    call ModeProducts(res, parts(:, 1), Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, transposeQ=.false.)
    call SetRealPart(f, res)

    if (size(parts, 2) > 1) then
      call ModeProducts(res, parts(:, 2), Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, transposeQ=.false.)
      call SetImagPart(f, res)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fem_SolveShiftedLaplacianReal(u, f, shift)
    use M_Grid

    complex(R64), intent(out), contiguous :: u(:)
    complex(R64), intent(in), contiguous  :: f(:)
    real(R64), intent(in)                 :: shift

    real(R64), allocatable :: parts(:, :)
    integer(I32) :: iPart

    call CheckFieldSize(size(u), "Grid_Fem_SolveShiftedLaplacian")
    call CheckFieldSize(size(f), "Grid_Fem_SolveShiftedLaplacian")
    call CheckShiftedSpectrum(shift, "Grid_Fem_SolveShiftedLaplacian")

    ! u = S3 (D + shift)^-1 S3^T W3 f, real and imaginary part separately
    call ExtractParts(parts, f, weightedQ=.true.)

    do iPart = 1, size(parts, 2)
      call ApplySpectralInverse(parts(:, iPart), shift)
    end do

    call SetRealPart(u, parts(:, 1))
    if (size(parts, 2) > 1) call SetImagPart(u, parts(:, 2))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fem_SolveShiftedLaplacianComplex(u, f, shift)
    use M_Grid

    complex(R64), intent(out), contiguous :: u(:)
    complex(R64), intent(in), contiguous  :: f(:)
    complex(R64), intent(in)              :: shift

    real(R64), allocatable :: part(:), cRe(:), cIm(:)

    ! A real shift takes the real path (half the work for a real f, and the
    ! check that D + shift stays positive)
    if (abs(aimag(shift)) <= 0.0_R64) then
      call Grid_Fem_SolveShiftedLaplacianReal(u, f, real(shift, kind=R64))
      return
    end if

    call CheckFieldSize(size(u), "Grid_Fem_SolveShiftedLaplacian")
    call CheckFieldSize(size(f), "Grid_Fem_SolveShiftedLaplacian")

    allocate (part(Grid_nPoints), cRe(Grid_nPoints), cIm(Grid_nPoints))

    ! c = S3^T W3 f, real transforms of the real and imaginary parts
    call ExtractPart(part, f, imagQ=.false., weightedQ=.true.)
    call ModeProducts(cRe, part, Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, transposeQ=.true.)

    if (HasImaginaryPart(f)) then
      call ExtractPart(part, f, imagQ=.true., weightedQ=.true.)
      call ModeProducts(cIm, part, Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, transposeQ=.true.)
    else
      cIm(:) = 0.0_R64
    end if

    ! c = c / (D + shift), where the two parts mix
    call DivideBySpectrumComplex(cRe, cIm, shift, Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz)

    ! u = S3 c
    call ModeProducts(part, cRe, Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, transposeQ=.false.)
    call SetRealPart(u, part)
    call ModeProducts(part, cIm, Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, transposeQ=.false.)
    call SetImagPart(u, part)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Mode products of the FDM eigenvectors with a real node array,
  !>   out = S3^T in   (transposeQ),   out = S3 in   (otherwise),
  !> as one DGEMM along x on the (nx, ny nz) matrix, nz DGEMMs along y on the
  !> (nx, ny) slices (OpenMP over z) and one DGEMM along z on the (nx ny, nz)
  !> matrix (in -> out -> tmp -> out).
  subroutine ModeProducts(out, in, nx, ny, nz, transposeQ)
    !$ use omp_lib, only: omp_in_parallel
    use M_Utils_BlasLib, only: DGEMM

    integer(I32), intent(in) :: nx
    integer(I32), intent(in) :: ny
    integer(I32), intent(in) :: nz
    real(R64), intent(out)   :: out(nx, ny, nz)
    real(R64), intent(in)    :: in(nx, ny, nz)
    logical, intent(in)      :: transposeQ

    real(R64), allocatable :: tmp(:, :, :)
    character(len=1) :: opX, opYZ
    integer(I32) :: iz

    ! S^T acts along x as op(Sx) = Sx^T from the left, along y and z as Sy and
    ! Sz from the right of the slices; S the other way round
    if (transposeQ) then
      opX = 'T'
      opYZ = 'N'
    else
      opX = 'N'
      opYZ = 'T'
    end if

    allocate (tmp(nx, ny, nz))

    call DGEMM(opX, 'N', nx, ny * nz, nx, 1.0_R64, axisOperators(1) % eigenvectors, nx, in, nx, 0.0_R64, out, nx)

    !$omp parallel do default(shared) private(iz) schedule(static) if (.not. omp_in_parallel())
    do iz = 1, nz
      call DGEMM('N', opYZ, nx, ny, ny, 1.0_R64, out(1, 1, iz), nx, axisOperators(2) % eigenvectors, ny, &
                 0.0_R64, tmp(1, 1, iz), nx)
    end do
    !$omp end parallel do

    call DGEMM('N', opYZ, nx * ny, nz, nz, 1.0_R64, tmp, nx * ny, axisOperators(3) % eigenvectors, nz, 0.0_R64, out, nx * ny)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> t = S3 (D + shift)^-1 S3^T t for a real node array (D + shift > 0).
  subroutine ApplySpectralInverse(t, shift)

    real(R64), intent(inout), contiguous :: t(:)
    real(R64), intent(in)                :: shift

    real(R64), allocatable :: c(:)

    allocate (c(size(t)))

    call ModeProducts(c, t, Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, transposeQ=.true.)
    call DivideBySpectrum(c, shift, Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz)
    call ModeProducts(t, c, Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, transposeQ=.false.)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> c(ix, iy, iz) = c(ix, iy, iz) / (lx(ix) + ly(iy) + lz(iz) + shift).
  subroutine DivideBySpectrum(c, shift, nx, ny, nz)
    !$ use omp_lib, only: omp_in_parallel

    integer(I32), intent(in) :: nx
    integer(I32), intent(in) :: ny
    integer(I32), intent(in) :: nz
    real(R64), intent(inout) :: c(nx, ny, nz)
    real(R64), intent(in)    :: shift

    real(R64) :: lyz
    integer(I32) :: ix, iy, iz

    !$omp parallel do default(shared) private(ix, iy, iz, lyz) schedule(static) if (.not. omp_in_parallel())
    do iz = 1, nz
      do iy = 1, ny
        lyz = Grid_Fem_axes(2) % eigenvalues(iy) + Grid_Fem_axes(3) % eigenvalues(iz) + shift
        do ix = 1, nx
          c(ix, iy, iz) = c(ix, iy, iz) / (Grid_Fem_axes(1) % eigenvalues(ix) + lyz)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> (cRe + i cIm)(ix, iy, iz) = (cRe + i cIm)(ix, iy, iz) / (lx(ix) + ly(iy) + lz(iz) + shift)
  !> for a complex shift (the coefficients are kept as two real arrays).
  subroutine DivideBySpectrumComplex(cRe, cIm, shift, nx, ny, nz)
    !$ use omp_lib, only: omp_in_parallel

    integer(I32), intent(in) :: nx
    integer(I32), intent(in) :: ny
    integer(I32), intent(in) :: nz
    real(R64), intent(inout) :: cRe(nx, ny, nz)
    real(R64), intent(inout) :: cIm(nx, ny, nz)
    complex(R64), intent(in) :: shift

    complex(R64) :: lyz, q
    integer(I32) :: ix, iy, iz

    !$omp parallel do default(shared) private(ix, iy, iz, lyz, q) schedule(static) if (.not. omp_in_parallel())
    do iz = 1, nz
      do iy = 1, ny
        lyz = Grid_Fem_axes(2) % eigenvalues(iy) + Grid_Fem_axes(3) % eigenvalues(iz) + shift
        do ix = 1, nx
          q = cmplx(cRe(ix, iy, iz), cIm(ix, iy, iz), kind=R64) / (Grid_Fem_axes(1) % eigenvalues(ix) + lyz)
          cRe(ix, iy, iz) = real(q, kind=R64)
          cIm(ix, iy, iz) = aimag(q)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> part = Re(f) or Im(f) (imagQ), times the 3D weights if weightedQ.
  subroutine ExtractPart(part, f, imagQ, weightedQ)
    !$ use omp_lib, only: omp_in_parallel

    real(R64), intent(out), contiguous   :: part(:)
    complex(R64), intent(in), contiguous :: f(:)
    logical, intent(in)                  :: imagQ
    logical, intent(in)                  :: weightedQ

    integer(I32) :: i

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(f)
      if (imagQ) then
        part(i) = aimag(f(i))
      else
        part(i) = real(f(i), kind=R64)
      end if
      if (weightedQ) part(i) = part(i) * Grid_Fem_weights(i)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> parts(:, 1) = Re(f) and, only if Im(f) is not identically zero,
  !> parts(:, 2) = Im(f) (size(parts, 2) = 1 or 2), times the 3D weights if
  !> weightedQ and times factor_ (default 1), in one pass over f.
  subroutine ExtractParts(parts, f, weightedQ, factor_)
    !$ use omp_lib, only: omp_in_parallel

    real(R64), allocatable, intent(out)  :: parts(:, :)
    complex(R64), intent(in), contiguous :: f(:)
    logical, intent(in)                  :: weightedQ
    real(R64), intent(in), optional      :: factor_

    real(R64) :: factor, scale
    integer(I32) :: i, nParts

    factor = 1.0_R64
    if (present(factor_)) factor = factor_

    nParts = 1
    if (HasImaginaryPart(f)) nParts = 2

    allocate (parts(size(f), nParts))

    !$omp parallel do default(shared) private(i, scale) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(f)
      scale = factor
      if (weightedQ) scale = scale * Grid_Fem_weights(i)
      parts(i, 1) = real(f(i), kind=R64) * scale
      if (nParts > 1) parts(i, 2) = aimag(f(i)) * scale
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> f = part + 0 i.
  subroutine SetRealPart(f, part)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(out), contiguous :: f(:)
    real(R64), intent(in), contiguous     :: part(:)

    integer(I32) :: i

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(f)
      f(i) = cmplx(part(i), 0.0_R64, kind=R64)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> f = Re(f) + part i.
  subroutine SetImagPart(f, part)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(inout), contiguous :: f(:)
    real(R64), intent(in), contiguous       :: part(:)

    integer(I32) :: i

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(f)
      f(i) = cmplx(real(f(i), kind=R64), part(i), kind=R64)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Whether the imaginary part of f is not identically zero (a NaN counts as
  !> nonzero, so that it propagates instead of being dropped).
  pure function HasImaginaryPart(f) result(res)

    complex(R64), intent(in), contiguous :: f(:)
    logical :: res

    integer(I32) :: i

    res = .false.
    do i = 1, size(f)
      if (.not. (abs(aimag(f(i))) <= 0.0_R64)) then
        res = .true.
        return
      end if
    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Stops unless a field argument has Grid_nPoints entries.
  subroutine CheckFieldSize(n, caller)
    use M_Grid

    integer(I32), intent(in)     :: n
    character(len=*), intent(in) :: caller

    if (n .ne. Grid_nPoints) error stop caller//": field size differs from Grid_nPoints"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Stops unless the shifted spectrum D + shift of -lap is positive, i.e.
  !> lx(1) + ly(1) + lz(1) + shift > 0 (smallest eigenvalues of the axes;
  !> NaN-safe), so that the spectral inverse (D + shift)^-1 exists.
  subroutine CheckShiftedSpectrum(shift, caller)

    real(R64), intent(in)        :: shift
    character(len=*), intent(in) :: caller

    if (.not. (Grid_Fem_axes(1) % eigenvalues(1) + Grid_Fem_axes(2) % eigenvalues(1) + &
               Grid_Fem_axes(3) % eigenvalues(1) + shift > 0.0_R64)) then
      error stop caller//": the shifted spectrum D + shift of -lap is not positive "// &
        "(a real shift must exceed -(smallest eigenvalue of -lap))"
    end if

  end subroutine

end submodule
