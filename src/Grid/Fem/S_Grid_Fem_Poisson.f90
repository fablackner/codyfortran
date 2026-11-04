! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Fem:S_Grid_Fem_Laplacian) S_Grid_Fem_Poisson
  !-----------------------------------------------------------------------------
  ! Galerkin Poisson solves of the FEM grid with free-space boundary values:
  ! Grid_Fem_SolvePoisson (multipole expansion of the load) and
  ! Grid_Fem_SolvePoissonPointCharges (exact values of point charges), and the
  ! tables of the multipole boundary values (built by Grid_Fem_SetupPoisson).
  ! A child of S_Grid_Fem_Laplacian (nested submodule): the boundary columns
  ! of its axis operator tables, the spectral inverse, the real and imaginary
  ! part helpers and the argument checks are reached by host association.
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Polynomial coefficients of the real regular solid harmonics (Racah
  !> normalization, M_Utils_SolidHarmonics) up to l = Grid_Fem_multipoleLmax,
  !> filled in Grid_Fem_SetupPoisson:
  !>   S_lm(x, y, z) = sum_(a+b+c=l) solidCoeffs(a, b, c, SolidHarmonics_LmIndex(l, m)) x^a y^b z^c,
  !> with the addition theorem 1/|r1 - r2| = sum_lm S_lm(r1) S_lm(r2)/|r2|^(2l+1)
  !> for |r1| < |r2| (S_10 = z, S_11 = x, S_1-1 = y, S_20 = (3z^2 - r^2)/2, ...).
  real(R64), allocatable :: solidCoeffs(:, :, :, :)

  !> Coordinates (3, nFaceNodes) of the box-surface nodes that couple to
  !> interior nodes (order of FillFaceNodes), filled in Grid_Fem_SetupPoisson.
  real(R64), allocatable :: faceNodes(:, :)

  !> Irregular solid harmonics of the face nodes r_B about the expansion
  !> center c, filled in Grid_Fem_SetupPoisson ((nFaceNodes, (lmax + 1)^2)):
  !>   faceHarmonics(iFace, SolidHarmonics_LmIndex(l, m)) = S_lm(r_B - c)/|r_B - c|^(2l+1),
  !> so that the multipole boundary values are V_B = faceHarmonics Q with the
  !> multipole moments Q_lm = sum_i load_i S_lm(r_i - c) (one DGEMV/DGEMM).
  real(R64), allocatable :: faceHarmonics(:, :)

  !> Number of divergence warnings of the multipole boundary values printed
  !> so far (at most one per run; atomic reads, written in a critical section).
  integer(I32) :: nMultipoleWarnings = 0

  !> Fraction of sum |load| outside the convergence sphere of the multipole
  !> boundary values above which Grid_Fem_SolvePoisson warns. The error of
  !> the boundary values is of the order of that fraction times the total
  !> load over the convergence radius R, so a tail of a bound density is
  !> harmless: 3.8e-7 (1s density, box +-10), 4.9e-7 (H2 RHF, box +-10) and
  !> 1.55e-5 (H2O with the diffuse STO start guess, box +-8) change the
  !> potential far below the discretization error. Loads that reach the
  !> faces (elongated boxes, ionized or off-center densities) have fractions
  !> of 1e-3 ... 1, with boundary-value errors up to O(1).
  real(R64), parameter :: MULTIPOLE_WARNING_FRACTION = 1.0e-4_R64

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fem_SetupPoisson
    use M_Utils_SolidHarmonics, only: SolidHarmonics_FillCoefficients

    integer(I32) :: iFace

    ! Multipole boundary values: solid-harmonic polynomials for the moments,
    ! face-node geometry, irregular solid harmonics at the face nodes and the
    ! convergence radius (distance of the expansion center to the nearest face node)
    call SolidHarmonics_FillCoefficients(solidCoeffs, Grid_Fem_multipoleLmax)
    call FillFaceNodes(faceNodes)
    call FillFaceHarmonics(faceHarmonics, faceNodes, Grid_Fem_multipoleCenter, Grid_Fem_multipoleLmax)

    Grid_Fem_multipoleRadius = huge(1.0_R64)
    do iFace = 1, size(faceNodes, 2)
      Grid_Fem_multipoleRadius = min(Grid_Fem_multipoleRadius, norm2(faceNodes(:, iFace) - Grid_Fem_multipoleCenter))
    end do

    ! Multipole lines of the mesh summary (after those of PrintMeshSummary)
    write (*, '(A, 3F10.4, A, F10.4)') "  multipole center", Grid_Fem_multipoleCenter, ", convergence radius", &
      Grid_Fem_multipoleRadius
    write (*, '(A, I0, A, I0, A, F0.1, A)') "  face harmonics table: ", size(faceHarmonics, 1), " face nodes x ", &
      size(faceHarmonics, 2), " (l, m) = ", 8.0_R64 * real(size(faceHarmonics, kind=I64), R64) / 1.0e6_R64, " MB"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fem_SolvePoisson(pot, load)
    use M_Utils_Constants, only: PI
    use M_Grid

    complex(R64), intent(out), contiguous :: pot(:)
    complex(R64), intent(in), contiguous  :: load(:)

    real(R64), allocatable :: parts(:, :), faceValues(:, :)
    integer(I32) :: iPart

    call CheckFieldSize(size(pot), "Grid_Fem_SolvePoisson")
    call CheckFieldSize(size(load), "Grid_Fem_SolvePoisson")
    call CheckShiftedSpectrum(0.0_R64, "Grid_Fem_SolvePoisson")

    ! Right-hand sides 4 pi Re(load) and, if nonzero, 4 pi Im(load)
    call ExtractParts(parts, load, weightedQ=.false., factor_=4.0_R64 * PI)
    call CheckMultipoleConvergence(load)

    ! Multipole boundary values of all parts (moments in one pass; the
    ! factor 4 pi of the parts is divided out)
    call FillMultipoleFaceValues(faceValues, parts, 1.0_R64 / (4.0_R64 * PI))

    ! K_II pot = 4 pi load - K_IB V_B for every part
    do iPart = 1, size(parts, 2)
      call SubtractBoundaryCoupling(parts(:, iPart), faceValues(:, iPart), Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz)
      call ApplySpectralInverse(parts(:, iPart), 0.0_R64)
    end do

    call SetRealPart(pot, parts(:, 1))
    if (size(parts, 2) > 1) call SetImagPart(pot, parts(:, 2))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Convergence diagnostic of the multipole boundary values: the expansion
  !> about c converges at the face nodes only if the whole load lies inside
  !> the sphere |r - c| < Grid_Fem_multipoleRadius (nearest face node). If
  !> more than MULTIPOLE_WARNING_FRACTION (1e-4) of sum |load| lies outside,
  !> a warning is printed, once per run (afterwards the check is skipped).
  subroutine CheckMultipoleConvergence(load)

    complex(R64), intent(in), contiguous :: load(:)

    real(R64) :: fraction
    integer(I32) :: nWarnings

    !$omp atomic read
    nWarnings = nMultipoleWarnings
    if (nWarnings > 0) return

    fraction = OutsideLoadFraction(load, Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz)
    if (.not. (fraction > MULTIPOLE_WARNING_FRACTION)) return

    !$omp critical (Grid_Fem_MultipoleWarning)
    if (nMultipoleWarnings .eq. 0) then
      write (*, '(A, ES9.2, A, F0.4, A, 3F9.4, A)') "  warning: grid.fem: a fraction ", fraction, &
        " of sum |load| of Grid_Fem_SolvePoisson lies outside the convergence sphere of the multipole boundary values (radius ", &
        Grid_Fem_multipoleRadius, " about", Grid_Fem_multipoleCenter, "): the potential may be inaccurate."
      write (*, '(A)') "  Center the box on the charge (and grid.fem.multipoleCenter) or enlarge the box "// &
        "(elongated boxes are affected most); this warning is printed once."
      !$omp atomic write
      nMultipoleWarnings = 1
    end if
    !$omp end critical (Grid_Fem_MultipoleWarning)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Fraction of sum_i |load_i| at the nodes outside the convergence sphere of
  !> the multipole expansion, |r_i - c| >= Grid_Fem_multipoleRadius (0 for a
  !> vanishing load).
  function OutsideLoadFraction(load, nx, ny, nz) result(fraction)
    !$ use omp_lib, only: omp_in_parallel

    integer(I32), intent(in) :: nx
    integer(I32), intent(in) :: ny
    integer(I32), intent(in) :: nz
    complex(R64), intent(in) :: load(nx, ny, nz)
    real(R64) :: fraction

    real(R64) :: dx2(nx), radius2, dz2, dyz2, a, sumAll, sumOutside
    integer(I32) :: ix, iy, iz

    radius2 = Grid_Fem_multipoleRadius**2
    dx2(:) = (Grid_Fem_axes(1) % points(:) - Grid_Fem_multipoleCenter(1))**2

    sumAll = 0.0_R64
    sumOutside = 0.0_R64

    !$omp parallel do default(shared) private(ix, iy, iz, dz2, dyz2, a) reduction(+: sumAll, sumOutside) &
    !$omp schedule(static) if (.not. omp_in_parallel())
    do iz = 1, nz
      dz2 = (Grid_Fem_axes(3) % points(iz) - Grid_Fem_multipoleCenter(3))**2
      do iy = 1, ny
        dyz2 = (Grid_Fem_axes(2) % points(iy) - Grid_Fem_multipoleCenter(2))**2 + dz2
        do ix = 1, nx
          a = abs(load(ix, iy, iz))
          sumAll = sumAll + a
          if (dx2(ix) + dyz2 >= radius2) sumOutside = sumOutside + a
        end do
      end do
    end do
    !$omp end parallel do

    fraction = 0.0_R64
    if (sumAll > 0.0_R64) fraction = sumOutside / sumAll

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fem_SolvePoissonPointCharges(pot, charges, positions)
    !$ use omp_lib, only: omp_in_parallel
    use M_Utils_Constants, only: PI
    use M_Grid

    complex(R64), intent(out), contiguous :: pot(:)
    real(R64), intent(in), contiguous     :: charges(:)
    real(R64), intent(in), contiguous     :: positions(:, :)

    real(R64), allocatable :: t(:), faceValues(:)
    integer(I32) :: iCharge, iAxis, iFace

    call CheckFieldSize(size(pot), "Grid_Fem_SolvePoissonPointCharges")
    call CheckShiftedSpectrum(0.0_R64, "Grid_Fem_SolvePoissonPointCharges")

    if (size(positions, 1) .ne. 3 .or. size(positions, 2) .ne. size(charges)) then
      error stop "Grid_Fem_SolvePoissonPointCharges: positions must have the shape (3, size(charges))"
    end if

    do iCharge = 1, size(charges)
      do iAxis = 1, 3
        if (.not. (positions(iAxis, iCharge) > Grid_Fem_axes(iAxis) % lowerBound .and. &
                   positions(iAxis, iCharge) < Grid_Fem_axes(iAxis) % upperBound)) then
          error stop "Grid_Fem_SolvePoissonPointCharges: every charge must lie strictly inside the box"
        end if
      end do
    end do

    ! Galerkin load 4 pi sum_I q_I phi_i(R_I)
    allocate (t(Grid_nPoints))
    t(:) = 0.0_R64
    do iCharge = 1, size(charges)
      call AddPointLoad(t, 4.0_R64 * PI * charges(iCharge), positions(:, iCharge), Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz)
    end do

    ! Exact box-surface values sum_I q_I/|r_B - R_I| at the face nodes
    allocate (faceValues(size(faceNodes, 2)))

    !$omp parallel do default(shared) private(iFace, iCharge) schedule(static) if (.not. omp_in_parallel())
    do iFace = 1, size(faceNodes, 2)
      faceValues(iFace) = 0.0_R64
      do iCharge = 1, size(charges)
        faceValues(iFace) = faceValues(iFace) + charges(iCharge) / norm2(faceNodes(:, iFace) - positions(:, iCharge))
      end do
    end do
    !$omp end parallel do

    call SubtractBoundaryCoupling(t, faceValues, Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz)
    call ApplySpectralInverse(t, 0.0_R64)
    call SetRealPart(pot, t)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Coordinates of the box-surface nodes that couple to interior nodes: face
  !> nodes whose two other coordinates are interior (edge and corner nodes do
  !> not couple through the lumped mass). Order: for the normal axis a = x, y,
  !> z the lower face, then the upper face, each with the nodes (kb, kc) of
  !> the tangential axes b < c (TangentialAxes), kb fastest: lower x face
  !> (iy, iz), upper x face, lower y face (ix, iz), upper y face, lower z face
  !> (ix, iy), upper z face; SubtractBoundaryCoupling reads this order.
  subroutine FillFaceNodes(nodes)

    real(R64), allocatable, intent(out) :: nodes(:, :)

    real(R64) :: node(3)
    integer(I32) :: n(3), a, b, c, iSide, kb, kc, iFace

    n = [Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz]
    allocate (nodes(3, 2 * (n(2) * n(3) + n(1) * n(3) + n(1) * n(2))))

    iFace = 0

    do a = 1, 3
      call TangentialAxes(b, c, a)
      do iSide = 1, 2
        node(a) = merge(Grid_Fem_axes(a) % lowerBound, Grid_Fem_axes(a) % upperBound, iSide .eq. 1)
        do kc = 1, n(c)
          node(c) = Grid_Fem_axes(c) % points(kc)
          do kb = 1, n(b)
            node(b) = Grid_Fem_axes(b) % points(kb)
            iFace = iFace + 1
            nodes(:, iFace) = node
          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> rhs = rhs - K_IB V_B for the box-surface values V_B (faceValues in the
  !> order of FillFaceNodes). A face node couples along its normal axis a only,
  !>   (K_IB V_B)(i) = Ka_full(i_a, box node) w_b(i_b) w_c(i_c) V_B(box node, i_b, i_c)
  !> (b < c the tangential axes, e.g. Kx_full(ix, box node) wy(iy) wz(iz) on
  !> the x faces); the boundary columns vanish outside the first/last element
  !> (p = nLocals - 1 interior nodes). Every rhs entry receives its terms in
  !> the order x, y, z faces, lower before upper.
  subroutine SubtractBoundaryCoupling(rhs, faceValues, nx, ny, nz)

    integer(I32), intent(in) :: nx
    integer(I32), intent(in) :: ny
    integer(I32), intent(in) :: nz
    real(R64), intent(inout) :: rhs(nx, ny, nz)
    real(R64), intent(in)    :: faceValues(:)

    real(R64) :: vLower, vUpper
    integer(I32) :: p, n(3), idx(3), a, b, c, kb, kc, k, nFace, offset, iFace

    p = Grid_Fem_nLocals - 1
    n = [nx, ny, nz]

    ! The faces normal to a: the lower one at offset, the upper one at offset + nFace
    offset = 0
    do a = 1, 3
      call TangentialAxes(b, c, a)
      nFace = n(b) * n(c)
      do kc = 1, n(c)
        idx(c) = kc
        do kb = 1, n(b)
          idx(b) = kb
          iFace = offset + kb + n(b) * (kc - 1)
          vLower = Grid_Fem_axes(b) % weights(kb) * Grid_Fem_axes(c) % weights(kc) * faceValues(iFace)
          vUpper = Grid_Fem_axes(b) % weights(kb) * Grid_Fem_axes(c) % weights(kc) * faceValues(iFace + nFace)
          do k = 1, min(p, n(a))
            idx(a) = k
            rhs(idx(1), idx(2), idx(3)) = rhs(idx(1), idx(2), idx(3)) - axisOperators(a) % lowerBoundaryColumn(k) * vLower
          end do
          do k = max(1, n(a) - p + 1), n(a)
            idx(a) = k
            rhs(idx(1), idx(2), idx(3)) = rhs(idx(1), idx(2), idx(3)) - axisOperators(a) % upperBoundaryColumn(k) * vUpper
          end do
        end do
      end do
      offset = offset + 2 * nFace
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> The tangential axes b < c of the box faces normal to axis a: (2, 3) for x,
  !> (1, 3) for y, (1, 2) for z.
  pure subroutine TangentialAxes(b, c, a)

    integer(I32), intent(out) :: b
    integer(I32), intent(out) :: c
    integer(I32), intent(in)  :: a

    b = merge(2, 1, a .eq. 1)
    c = merge(2, 3, a .eq. 3)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Multipole boundary values of the real loads parts(:, k), k = 1..nParts, at the face nodes,
  !>   faceValues(:, k) = scale faceHarmonics Q(:, k),   Q_lm = sum_i parts(i, k) S_lm(r_i - c),
  !> one DGEMV (one part) or DGEMM. The moments of all parts come from one pass
  !> over parts: the Cartesian moments M(a, b, c) = sum_i parts_i x^a y^b z^c
  !> (r_i - c = (x, y, z); one DGEMM along x for all parts, then contractions
  !> with the y and z power tables) and Q_lm = sum_abc s_lm(a, b, c) M(a, b, c)
  !> (s_lm = solidCoeffs, one DGEMM).
  subroutine FillMultipoleFaceValues(faceValues, parts, scale)
    use M_Utils_BlasLib, only: DGEMM, DGEMV

    real(R64), allocatable, intent(out) :: faceValues(:, :)
    real(R64), intent(in), contiguous   :: parts(:, :)
    real(R64), intent(in)               :: scale

    real(R64), allocatable :: powersX(:, :), powersY(:, :), powersZ(:, :)
    real(R64), allocatable :: momentsX(:, :, :, :), momentsXY(:, :, :, :), moments(:, :, :, :), multipoles(:, :)
    integer(I32) :: lmax, nx, ny, nz, nParts, nFaces, nLm, nCartesian, iPart, iy, iz, b, c

    lmax = Grid_Fem_multipoleLmax
    nx = Grid_Fem_nx
    ny = Grid_Fem_ny
    nz = Grid_Fem_nz
    nParts = size(parts, 2)
    nFaces = size(faceHarmonics, 1)
    nLm = (lmax + 1)**2
    nCartesian = (lmax + 1)**3

    call FillPowerTable(powersX, Grid_Fem_axes(1) % points - Grid_Fem_multipoleCenter(1), lmax)
    call FillPowerTable(powersY, Grid_Fem_axes(2) % points - Grid_Fem_multipoleCenter(2), lmax)
    call FillPowerTable(powersZ, Grid_Fem_axes(3) % points - Grid_Fem_multipoleCenter(3), lmax)

    ! momentsX(a, iy, iz, k) = sum_ix x^a parts(ix, iy, iz, k), all parts at once
    allocate (momentsX(0:lmax, ny, nz, nParts))
    call DGEMM('N', 'N', lmax + 1, ny * nz * nParts, nx, 1.0_R64, powersX, lmax + 1, parts, nx, 0.0_R64, momentsX, &
               lmax + 1)

    ! momentsXY(a, b, iz, k) = sum_iy momentsX(a, iy, iz, k) y^b
    allocate (momentsXY(0:lmax, 0:lmax, nz, nParts))
    momentsXY = 0.0_R64
    do iPart = 1, nParts
      do iz = 1, nz
        do iy = 1, ny
          do b = 0, lmax
            momentsXY(:, b, iz, iPart) = momentsXY(:, b, iz, iPart) + momentsX(:, iy, iz, iPart) * powersY(b, iy)
          end do
        end do
      end do
    end do

    ! moments(a, b, c, k) = sum_iz momentsXY(a, b, iz, k) z^c
    allocate (moments(0:lmax, 0:lmax, 0:lmax, nParts))
    moments = 0.0_R64
    do iPart = 1, nParts
      do iz = 1, nz
        do c = 0, lmax
          moments(:, :, c, iPart) = moments(:, :, c, iPart) + momentsXY(:, :, iz, iPart) * powersZ(c, iz)
        end do
      end do
    end do

    ! Multipole moments Q(lm, k) = sum_abc solidCoeffs(a, b, c, lm) moments(a, b, c, k)
    allocate (multipoles(nLm, nParts))
    call DGEMM('T', 'N', nLm, nParts, nCartesian, 1.0_R64, solidCoeffs, nCartesian, moments, nCartesian, 0.0_R64, &
               multipoles, nLm)

    ! V_B = scale faceHarmonics Q
    allocate (faceValues(nFaces, nParts))
    if (nParts .eq. 1) then
      call DGEMV('N', nFaces, nLm, scale, faceHarmonics, nFaces, multipoles, 1, 0.0_R64, faceValues, 1)
    else
      call DGEMM('N', 'N', nFaces, nParts, nLm, scale, faceHarmonics, nFaces, multipoles, nLm, 0.0_R64, faceValues, nFaces)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Powers table(a, i) = coords(i)^a, a = 0..lmax.
  subroutine FillPowerTable(table, coords, lmax)

    real(R64), allocatable, intent(out) :: table(:, :)
    real(R64), intent(in)               :: coords(:)
    integer(I32), intent(in)            :: lmax

    integer(I32) :: a

    allocate (table(0:lmax, size(coords)))
    table(0, :) = 1.0_R64
    do a = 1, lmax
      table(a, :) = table(a - 1, :) * coords(:)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Irregular solid harmonics of the face nodes about the expansion center,
  !>   table(iFace, SolidHarmonics_LmIndex(l, m)) = S_lm(r_B - center)/|r_B - center|^(2l+1),
  !> 0 <= l <= lmax, r_B = nodes(:, iFace) (/= center): the boundary values of
  !> the multipole expansion per unit moment Q_lm (addition theorem).
  subroutine FillFaceHarmonics(table, nodes, center, lmax)
    !$ use omp_lib, only: omp_in_parallel
    use M_Utils_SolidHarmonics, only: SolidHarmonics_FillValues, SolidHarmonics_LmIndex

    real(R64), allocatable, intent(out) :: table(:, :)
    real(R64), intent(in)               :: nodes(:, :)
    real(R64), intent(in)               :: center(3)
    integer(I32), intent(in)            :: lmax

    real(R64) :: values((lmax + 1)**2), r(3), rInv, radial
    integer(I32) :: iFace, l, m, lm

    allocate (table(size(nodes, 2), (lmax + 1)**2))

    !$omp parallel do default(shared) private(iFace, l, m, lm, values, r, rInv, radial) schedule(static) &
    !$omp if (.not. omp_in_parallel())
    do iFace = 1, size(nodes, 2)
      r(:) = nodes(:, iFace) - center(:)
      call SolidHarmonics_FillValues(values, r, lmax)

      rInv = 1.0_R64 / norm2(r)
      radial = rInv
      do l = 0, lmax
        do m = -l, l
          lm = SolidHarmonics_LmIndex(l, m)
          table(iFace, lm) = values(lm) * radial
        end do
        radial = radial * rInv * rInv
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> t = t + q phi_i(r) for the tensor Lagrange basis functions phi_i of the
  !> interior nodes; only the (p + 1)^3 nodes of the element containing r
  !> contribute (q delta_iJ if r is the node J).
  subroutine AddPointLoad(t, q, r, nx, ny, nz)

    integer(I32), intent(in) :: nx
    integer(I32), intent(in) :: ny
    integer(I32), intent(in) :: nz
    real(R64), intent(inout) :: t(nx, ny, nz)
    real(R64), intent(in)    :: q
    real(R64), intent(in)    :: r(3)

    real(R64), allocatable :: values(:, :)
    integer(I32) :: p, iAxis, first(3), kx, ky, kz, ix, iy, iz

    p = Grid_Fem_nLocals - 1

    allocate (values(0:p, 3))
    do iAxis = 1, 3
      call FillLocalBasisValues(values(:, iAxis), first(iAxis), Grid_Fem_axes(iAxis), r(iAxis), p)
    end do

    ! Box nodes (index 0 or n + 1) are Dirichlet nodes without a basis function
    do kz = 0, p
      iz = first(3) + kz
      if (iz < 1 .or. iz > nz) cycle
      do ky = 0, p
        iy = first(2) + ky
        if (iy < 1 .or. iy > ny) cycle
        do kx = 0, p
          ix = first(1) + kx
          if (ix < 1 .or. ix > nx) cycle
          t(ix, iy, iz) = t(ix, iy, iz) + q * values(kx, 1) * values(ky, 2) * values(kz, 3)
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Lagrange basis values L_k(x), k = 0..p, of the element containing x (the
  !> left one on a shared element boundary, where the C0 basis is single
  !> valued) and the interior index first of its local node 0; the local node
  !> k has the interior index first + k (0 and nPoints + 1: box nodes).
  subroutine FillLocalBasisValues(values, first, axis, x, p)

    integer(I32), intent(in)              :: p
    real(R64), intent(out)                :: values(0:p)
    integer(I32), intent(out)             :: first
    type(T_Grid_Fem_Axis), intent(in)     :: axis
    real(R64), intent(in)                 :: x

    real(R64) :: nodes(0:p)
    integer(I32) :: iE, k, j

    iE = 1
    do while (iE < axis % nElements .and. x > axis % bounds(iE))
      iE = iE + 1
    end do
    first = (iE - 1) * p

    do k = 0, p
      if (first + k .eq. 0) then
        nodes(k) = axis % lowerBound
      else if (first + k .eq. axis % nPoints + 1) then
        nodes(k) = axis % upperBound
      else
        nodes(k) = axis % points(first + k)
      end if
    end do

    do k = 0, p
      values(k) = 1.0_R64
      do j = 0, p
        if (j .ne. k) values(k) = values(k) * (x - nodes(j)) / (nodes(k) - nodes(j))
      end do
    end do

  end subroutine

end submodule
