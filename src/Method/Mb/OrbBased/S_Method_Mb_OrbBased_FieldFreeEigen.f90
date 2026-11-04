! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased) S_Method_Mb_OrbBased_FieldFreeEigen
  !-----------------------------------------------------------------------------
  ! Exact propagation with the field-free one-body Hamiltonian h0 = T + V_ext
  ! on Ylm grids (method.mb.orbBased.exactFieldFreeQ): the *Eigen procedures,
  ! which Method_Mb_OrbBased_Fabricate binds to the exact-part pointers
  ! (Method_Mb_OrbBased_SetupFieldFree, _PropagateFieldFree,
  ! _ApplyFieldFreePhi, _ProjectBelowFieldFreeCut) for grid.ylm. Other grids:
  ! S_Method_Mb_OrbBased_FieldFreeCayley (with the frozen mean-field reference
  ! of S_Method_Mb_OrbBased_FieldFreeMeanField).
  !
  ! For a central potential h0 is block diagonal in (l,m) and independent of
  ! m. In the metric-symmetrized radial basis (w^{1/2} f, w the radial metric
  ! weights including r^2) every block is real symmetric, hs_l = U_l E_l U_l^T,
  ! so that
  !   exp(-i h0 dt) f_lm = w^{-1/2} U_l exp(-i E_l dt) U_l^T w^{1/2} f_lm.
  ! The blocks are built once by applying the kinetic and potential operators
  ! to radial unit vectors placed in all channels at once. Exact propagation
  ! removes the stiffness of h0 (centrifugal barrier and kinetic energy near
  ! the origin) from the Krylov integrator of the linear split step.
  !
  ! Spectral cut (method.mb.orbBased.exactFieldFreeCut = Ec): only the nKeep_l
  ! eigenvectors with E < Ec are kept. The highest FEDVR eigenvalues (up to
  ! ~1e5 a.u., from the dense Lobatto points at the element boundaries) are
  ! unphysical; restricting the dynamics to P = sum_{E < Ec} |u><u| (exact
  ! part on the kept vectors, P applied to the orbital derivatives) removes
  ! them as a source of stiffness.
  !
  ! With declared m (method.mb.orbBased.mOfOrbs) every orbital fills only the
  ! column of its m in the block of l (blockColumns, built once with the
  ! blocks); the columns of other m stay exactly zero and are not transformed.
  !-----------------------------------------------------------------------------

  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Eigenvectors U_l of the symmetrized blocks, (nRadial, nRadial, 0:lmax).
  real(R64), allocatable :: blockEvecs(:, :, :)

  !> Eigenvalues E_l of the symmetrized blocks, (nRadial, 0:lmax).
  real(R64), allocatable :: blockEvals(:, :)

  !> Square roots of the radial metric weights.
  real(R64), allocatable :: sqrtWeights(:)

  !> Number of kept eigenvectors per l (all without spectral cut), (0:lmax).
  integer(I32), allocatable :: nKeep(:)

  !> (orbital, m) columns of the block of l that are transformed: the column
  !> of the declared m of every orbital, or all columns without declared m;
  !> nColMax is the largest count over l.
  type :: T_BlockColumns
    integer(I32), allocatable :: orb(:)
    integer(I32), allocatable :: m(:)
  end type
  type(T_BlockColumns), allocatable :: blockColumns(:)
  integer(I32) :: nColMax = 0

  !> GPU variant (method.mb.orbBased.exactFieldFreeGpuQ): eigenvectors on the
  !> device, device work arrays and the host staging buffer (2 nColMax
  !> columns per l).
  logical :: gpuQ = .false.
  type(c_ptr) :: deviceEvecs = c_null_ptr
  type(c_ptr) :: deviceX = c_null_ptr
  type(c_ptr) :: deviceY = c_null_ptr
  real(R64), allocatable, target :: hostBlocks(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_PropagateFieldFreeEigen(orbs, dt, imagTimeQ)
    use M_Utils_BlasLib, only: DGEMM
    use M_Grid_Ylm

    complex(R64), intent(inout), contiguous :: orbs(:, :)
    real(R64), intent(in) :: dt
    logical, intent(in) :: imagTimeQ

    integer(I32) :: nRad, l, nCol, nK
    real(R64), allocatable :: x(:, :), y(:, :)

    if (gpuQ) then
      call PropagateFieldFreeGpu(orbs, dt, imagTimeQ)
      return
    end if

    nRad = Grid_Ylm_nRadial

    !$omp parallel do default(shared) schedule(dynamic) private(l, nCol, nK, x, y)
    do l = 0, Grid_Ylm_lmax
      nCol = size(blockColumns(l) % orb)
      if (nCol .eq. 0) cycle
      nK = nKeep(l)
      allocate (x(nRad, 2 * nCol), y(nK, 2 * nCol))

      ! Eigenbasis coefficients U^T x, phases, back transformation U y (kept
      ! eigenvectors only, which also projects onto them)
      call GatherBlock(x, orbs, l)
      call DGEMM('T', 'N', nK, 2 * nCol, nRad, 1.0_R64, blockEvecs(:, :, l), nRad, x, nRad, 0.0_R64, y, nK)
      call ApplyPhases(y, l, dt, imagTimeQ)
      call DGEMM('N', 'N', nRad, 2 * nCol, nK, 1.0_R64, blockEvecs(:, :, l), nRad, y, nK, 0.0_R64, x, nRad)
      call ScatterBlock(orbs, x, l)

      deallocate (x, y)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_ApplyFieldFreePhiEigen(orbs, orbsList, dt)
    !---------------------------------------------------------------------------
    ! Per block: every slice is transformed to the eigenbasis, weighted with
    ! phi_m(-i E dt) and summed; one back transformation (CPU only, the host
    ! keeps the eigenvectors also in the GPU variant).
    !---------------------------------------------------------------------------
    use M_Utils_BlasLib, only: DGEMM
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: orbsList(:, :, 0:)
    real(R64), intent(in) :: dt

    integer(I32) :: nRad, nPhi, l, nCol, nK, m, k
    real(R64), allocatable :: x(:, :), y(:, :)
    complex(R64), allocatable :: phi(:, :), acc(:, :)

    nRad = Grid_Ylm_nRadial
    nPhi = size(orbsList, 3)
    orbs(:, :) = (0.0_R64, 0.0_R64)

    !$omp parallel do default(shared) schedule(dynamic) private(l, nCol, nK, m, k, x, y, phi, acc)
    do l = 0, Grid_Ylm_lmax
      nCol = size(blockColumns(l) % orb)
      if (nCol .eq. 0) cycle
      nK = nKeep(l)
      allocate (x(nRad, 2 * nCol), y(nK, 2 * nCol), phi(nK, 0:nPhi - 1), acc(nK, nCol))

      do k = 1, nK
        phi(k, :) = PhiFunctions(cmplx(0.0_R64, -blockEvals(k, l) * dt, kind=R64), nPhi - 1)
      end do

      acc(:, :) = (0.0_R64, 0.0_R64)
      do m = 0, nPhi - 1
        call GatherBlock(x, orbsList(:, :, m), l)
        call DGEMM('T', 'N', nK, 2 * nCol, nRad, 1.0_R64, blockEvecs(:, :, l), nRad, x, nRad, 0.0_R64, y, nK)
        acc(:, :) = acc(:, :) + spread(phi(:, m), 2, nCol) * cmplx(y(:, 1:nCol), y(:, nCol + 1:), kind=R64)
      end do

      y(:, 1:nCol) = real(acc, kind=R64)
      y(:, nCol + 1:) = aimag(acc)
      call DGEMM('N', 'N', nRad, 2 * nCol, nK, 1.0_R64, blockEvecs(:, :, l), nRad, y, nK, 0.0_R64, x, nRad)
      call ScatterBlock(orbs, x, l)

      deallocate (x, y, phi, acc)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> phi_0(z) ... phi_K(z): Taylor series for |z| < 1 (the recurrence
  !> phi_m = (phi_{m-1} - 1/(m-1)!) / z cancels there), the recurrence otherwise.
  pure function PhiFunctions(z, nMax) result(res)
    complex(R64), intent(in) :: z
    integer(I32), intent(in) :: nMax
    complex(R64) :: res(0:nMax)

    integer(I32), parameter :: nTaylor = 30
    integer(I32) :: m, j
    complex(R64) :: term
    real(R64) :: invFactorial

    if (abs(z) < 1.0_R64) then
      do m = 0, nMax
        ! sum_j z^j / (j + m)!
        invFactorial = 1.0_R64 / gamma(real(m + 1, kind=R64))
        term = cmplx(invFactorial, 0.0_R64, kind=R64)
        res(m) = term
        do j = 1, nTaylor
          term = term * z / real(j + m, kind=R64)
          res(m) = res(m) + term
        end do
      end do
    else
      res(0) = exp(z)
      invFactorial = 1.0_R64
      do m = 1, nMax
        res(m) = (res(m - 1) - invFactorial) / z
        invFactorial = invFactorial / real(m, kind=R64)
      end do
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine PropagateFieldFreeGpu(orbs, dt, imagTimeQ)
    !---------------------------------------------------------------------------
    ! GPU variant (method.mb.orbBased.exactFieldFreeGpuQ): the eigenvectors
    ! stay on the device; per call all blocks are gathered on the host, sent
    ! to the device once, transformed by one cuBLAS GEMM per l, and the phases
    ! are applied on the host between the two transformations.
    !---------------------------------------------------------------------------
    use M_Utils_CublasLib
    use M_Grid_Ylm

    complex(R64), intent(inout), contiguous :: orbs(:, :)
    real(R64), intent(in) :: dt
    logical, intent(in) :: imagTimeQ

    integer(I32) :: nRad, lmax, l, nCol
    real(R64), contiguous, pointer :: blocks(:, :, :)

    nRad = Grid_Ylm_nRadial
    lmax = Grid_Ylm_lmax

    blocks(1:nRad, 1:2 * nColMax, 0:lmax) => hostBlocks

    !$omp parallel do default(shared) schedule(dynamic) private(l, nCol)
    do l = 0, lmax
      nCol = size(blockColumns(l) % orb)
      call GatherBlock(blocks(:, 1:2 * nCol, l), orbs, l)
    end do
    !$omp end parallel do

    call TransformBlocksGpu('T', deviceX, deviceY)

    !$omp parallel do default(shared) schedule(dynamic) private(l, nCol)
    do l = 0, lmax
      nCol = size(blockColumns(l) % orb)
      call ApplyPhases(blocks(:, 1:2 * nCol, l), l, dt, imagTimeQ)
    end do
    !$omp end parallel do

    call TransformBlocksGpu('N', deviceY, deviceX)

    !$omp parallel do default(shared) schedule(dynamic) private(l, nCol)
    do l = 0, lmax
      nCol = size(blockColumns(l) % orb)
      call ScatterBlock(orbs, blocks(:, 1:2 * nCol, l), l)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> One cuBLAS transformation of all blocks: uploads hostBlocks to deviceIn,
  !> computes op(U_l) x_l per l into deviceOut, and downloads deviceOut to
  !> hostBlocks.
  subroutine TransformBlocksGpu(transA, deviceIn, deviceOut)
    use M_Utils_CublasLib
    use M_Grid_Ylm

    character(len=1), intent(in) :: transA
    type(c_ptr), intent(in) :: deviceIn
    type(c_ptr), intent(in) :: deviceOut

    integer(I32) :: nRad, lmax, l, nCol
    integer(I64) :: blockSize, offsetU, offsetX

    nRad = Grid_Ylm_nRadial
    lmax = Grid_Ylm_lmax
    blockSize = int(nRad, kind=I64) * 2 * nColMax

    call CublasLib_CopyToDevice(deviceIn, hostBlocks, 0_I64)
    do l = 0, lmax
      nCol = size(blockColumns(l) % orb)
      offsetU = int(nRad, kind=I64) * nRad * l
      offsetX = blockSize * l
      ! kept eigenvectors only: U_K^T x has nKeep rows (the rows below are
      ! not used by the back transformation U_K y)
      if (transA .eq. 'T') then
        call CublasLib_Dgemm('T', 'N', nKeep(l), 2 * nCol, nRad, 1.0_R64, deviceEvecs, offsetU, nRad, &
                             deviceIn, offsetX, nRad, 0.0_R64, deviceOut, offsetX, nRad)
      else
        call CublasLib_Dgemm('N', 'N', nRad, 2 * nCol, nKeep(l), 1.0_R64, deviceEvecs, offsetU, nRad, &
                             deviceIn, offsetX, nRad, 0.0_R64, deviceOut, offsetX, nRad)
      end if
    end do
    call CublasLib_CopyToHost(hostBlocks, deviceOut, 0_I64)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Gathers w^{1/2} f_lm of the (orbital, m) columns of the block of l
  !> (blockColumns): real parts in columns 1:nCol, imaginary parts in
  !> nCol+1:2 nCol (one GEMM per transformation).
  subroutine GatherBlock(x, orbs, l)
    use M_Grid_Ylm

    real(R64), intent(out), contiguous :: x(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    integer(I32), intent(in) :: l

    integer(I32) :: jCol, nCol
    complex(R64), allocatable :: fLm(:)

    allocate (fLm(Grid_Ylm_nRadial))
    associate (cols => blockColumns(l))
      nCol = size(cols % orb)
      do jCol = 1, nCol
        call Grid_Ylm_GetLmComponent(fLm, l, cols % m(jCol), orbs(:, cols % orb(jCol)))
        x(:, jCol) = sqrtWeights(:) * real(fLm(:), kind=R64)
        x(:, nCol + jCol) = sqrtWeights(:) * aimag(fLm(:))
      end do
    end associate

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Multiplies the eigenbasis coefficients of the block of l by
  !> exp(-i E dt) (or exp(-E dt) in imaginary time).
  subroutine ApplyPhases(y, l, dt, imagTimeQ)
    real(R64), intent(inout), contiguous :: y(:, :)
    integer(I32), intent(in) :: l
    real(R64), intent(in) :: dt
    logical, intent(in) :: imagTimeQ

    integer(I32) :: iCol, nCol, k
    complex(R64), allocatable :: phase(:)
    complex(R64) :: c

    allocate (phase(size(y, 1)))
    associate (evals => blockEvals(1:size(y, 1), l))
      if (imagTimeQ) then
        phase(:) = cmplx(exp(-evals * dt), 0.0_R64, kind=R64)
      else
        phase(:) = cmplx(cos(evals * dt), -sin(evals * dt), kind=R64)
      end if
    end associate

    nCol = size(y, 2) / 2
    do iCol = 1, nCol
      do k = 1, size(y, 1)
        c = cmplx(y(k, iCol), y(k, nCol + iCol), kind=R64) * phase(k)
        y(k, iCol) = real(c, kind=R64)
        y(k, nCol + iCol) = aimag(c)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_ProjectBelowFieldFreeCutEigen(orbs)
    !---------------------------------------------------------------------------
    ! orbs <- P orbs with P = 1 - U_H U_H^T per block in the symmetrized basis,
    ! U_H the eigenvectors above the cut (few compared to the kept ones).
    !---------------------------------------------------------------------------
    use M_Utils_BlasLib, only: DGEMM
    use M_Grid_Ylm

    complex(R64), intent(inout), contiguous :: orbs(:, :)

    integer(I32) :: nRad, l, nCol, nK, nH
    real(R64), allocatable :: x(:, :), y(:, :)

    nRad = Grid_Ylm_nRadial

    !$omp parallel do default(shared) schedule(dynamic) private(l, nCol, nK, nH, x, y)
    do l = 0, Grid_Ylm_lmax
      nK = nKeep(l)
      nH = nRad - nK
      if (nH .eq. 0) cycle
      nCol = size(blockColumns(l) % orb)
      if (nCol .eq. 0) cycle
      allocate (x(nRad, 2 * nCol), y(nH, 2 * nCol))

      call GatherBlock(x, orbs, l)
      call DGEMM('T', 'N', nH, 2 * nCol, nRad, 1.0_R64, blockEvecs(:, nK + 1:, l), nRad, x, nRad, 0.0_R64, y, nH)
      call DGEMM('N', 'N', nRad, 2 * nCol, nH, -1.0_R64, blockEvecs(:, nK + 1:, l), nRad, y, nH, 1.0_R64, x, nRad)
      call ScatterBlock(orbs, x, l)

      deallocate (x, y)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Scatters w^{-1/2} x of the columns of the block of l back into the
  !> channels (the other columns are left unchanged).
  subroutine ScatterBlock(orbs, x, l)
    use M_Grid_Ylm

    complex(R64), intent(inout), contiguous :: orbs(:, :)
    real(R64), intent(in), contiguous :: x(:, :)
    integer(I32), intent(in) :: l

    integer(I32) :: jCol, nCol
    complex(R64), allocatable :: fLm(:)

    allocate (fLm(Grid_Ylm_nRadial))
    associate (cols => blockColumns(l))
      nCol = size(cols % orb)
      do jCol = 1, nCol
        fLm(:) = cmplx(x(:, jCol), x(:, nCol + jCol), kind=R64) / sqrtWeights(:)
        call Grid_Ylm_SetLmComponent(orbs(:, cols % orb(jCol)), l, cols % m(jCol), fLm)
      end do
    end associate

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_SetupFieldFreeEigen
    !---------------------------------------------------------------------------
    ! Builds and diagonalizes the symmetrized radial blocks of h0 for every l.
    !
    ! Probe i is the radial unit vector e_i placed in every (l,m) channel; since
    ! h0 does not couple channels, channel (l,m) of h0 e_i is column i of the
    ! block of l. The m != 0 channels must reproduce the m = 0 block, and the
    ! symmetrized blocks must be real symmetric (Hermitian metric).
    !---------------------------------------------------------------------------
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_LapackLib
    use M_Grid
    use M_Grid_Ylm
    use M_SysKinetic
    use M_SysPotential
    use M_SysPotential_Ylm
    use M_Utils_CublasLib
    use M_Orbs

    real(R64), parameter :: tolerance = 1.0e-10_R64

    integer(I32) :: nRad, nG, nOS, lmax, mmax, i, l, m, nFound, iOrb, mCap
    integer(I64) :: blockSize
    logical, allocatable :: keepQ(:, :)
    logical :: potentialQ
    real(R64) :: hScale, devM, devImag, devSym
    complex(R64), allocatable :: probe(:), hProbe(:), tmp(:), hLm(:), unitRad(:), externalPotential(:)
    real(R64), allocatable :: hs(:, :), evals(:), evecs(:, :)

    call Say_Setup("method.mb.orbBased.exactFieldFree")

    !------------------------------------
    ! applicability
    !------------------------------------

    if (.not. Json_GetExistence("grid.ylm")) &
      error stop "method.mb.orbBased.exactFieldFreeQ requires a Ylm grid"
    if (Grid_Ylm_cProductQ) &
      error stop "method.mb.orbBased.exactFieldFreeQ requires a Hermitian (non-ECS) Ylm grid"
    if (.not. (SysKinetic_timeIndependentQ .and. SysKinetic_bodyTypeIndependentQ)) &
      error stop "method.mb.orbBased.exactFieldFreeQ requires a time- and body-type-independent kinetic operator"

    potentialQ = Json_GetExistence("sysPotential")
    if (potentialQ) then
      if (.not. (SysPotential_timeIndependentQ .and. SysPotential_bodyTypeIndependentQ)) &
        error stop "method.mb.orbBased.exactFieldFreeQ requires a time- and body-type-independent potential"
      if (.not. Json_GetExistence("sysPotential.ylm") .or. SysPotential_Ylm_lmax .ne. 0) &
        error stop "method.mb.orbBased.exactFieldFreeQ requires a central (l = 0) external potential"
    end if

    !------------------------------------
    ! probe h0 with radial unit vectors
    !------------------------------------

    nRad = Grid_Ylm_nRadial
    nG = Grid_nPoints
    lmax = Grid_Ylm_lmax
    mmax = Grid_Ylm_mmax

    allocate (sqrtWeights(nRad))
    sqrtWeights(:) = sqrt(Grid_Ylm_radialWeights(:))

    allocate (probe(nG), hProbe(nG), tmp(nG), hLm(nRad), unitRad(nRad))
    allocate (blockEvecs(nRad, nRad, 0:lmax), blockEvals(nRad, 0:lmax))

    if (potentialQ) call SysPotential_FillExternalPotential(externalPotential, 0.0_R64)

    devM = 0.0_R64
    devImag = 0.0_R64

    do i = 1, nRad
      probe(:) = (0.0_R64, 0.0_R64)
      unitRad(:) = (0.0_R64, 0.0_R64)
      unitRad(i) = (1.0_R64, 0.0_R64)
      do l = 0, lmax
        do m = -min(l, mmax), min(l, mmax)
          call Grid_Ylm_SetLmComponent(probe, l, m, unitRad)
        end do
      end do

      hProbe(:) = (0.0_R64, 0.0_R64)
      call SysKinetic_MultiplyWithKineticOp(hProbe, probe, 0.0_R64, bt_=1)
      if (potentialQ) then
        call SysPotential_MultiplyWithExternalPotential(tmp, externalPotential, probe)
        hProbe(:) = hProbe(:) + tmp(:)
      end if

      ! Column i of the block of l (from m = 0) in metric-symmetrized form
      do l = 0, lmax
        call Grid_Ylm_GetLmComponent(hLm, l, 0, hProbe)
        devImag = max(devImag, maxval(abs(aimag(hLm))))
        blockEvecs(:, i, l) = sqrtWeights(:) * real(hLm(:), kind=R64) / sqrtWeights(i)

        do m = -min(l, mmax), min(l, mmax)
          if (m .eq. 0) cycle
          call Grid_Ylm_GetLmComponent(tmp(1:nRad), l, m, hProbe)
          devM = max(devM, maxval(abs(tmp(1:nRad) - hLm(:))))
        end do
      end do
    end do

    !------------------------------------
    ! checks and diagonalization
    !------------------------------------

    hScale = maxval(abs(blockEvecs))
    devImag = devImag / hScale
    devM = devM / hScale

    devSym = 0.0_R64
    do l = 0, lmax
      devSym = max(devSym, maxval(abs(blockEvecs(:, :, l) - transpose(blockEvecs(:, :, l)))) / hScale)
    end do

    write (*, '(A,ES10.2,A,ES10.2,A,ES10.2)') " exact field-free blocks: asymmetry", devSym, &
      "  imaginary part", devImag, "  m dependence", devM

    if (devImag > tolerance) error stop "method.mb.orbBased.exactFieldFreeQ: h0 is not real on this grid"
    if (devM > tolerance) error stop "method.mb.orbBased.exactFieldFreeQ: h0 depends on m"
    if (devSym > tolerance) error stop "method.mb.orbBased.exactFieldFreeQ: h0 is not Hermitian in the grid metric"

    ! MRRR (DSYEVR): the QR iteration of DSYEV took most of the setup time
    !$omp parallel do default(shared) schedule(dynamic) private(l, hs, evals, evecs, nFound)
    do l = 0, lmax
      hs = 0.5_R64 * (blockEvecs(:, :, l) + transpose(blockEvecs(:, :, l)))
      call LapackLib_DiagonalizeSym(evals, evecs, nFound, hs, .true.)
      if (nFound .ne. nRad) error stop "method.mb.orbBased.exactFieldFreeQ: incomplete diagonalization of h0"
      blockEvals(:, l) = evals(:)
      blockEvecs(:, :, l) = evecs(:, :)
    end do
    !$omp end parallel do

    write (*, '(A,I0,A,ES12.4,A,ES12.4)') " exact field-free blocks diagonalized for l <= ", lmax, &
      ": eigenvalues from", minval(blockEvals), " to", maxval(blockEvals)

    ! Spectral cut: the eigenvalues are ascending, keep the first nKeep_l
    allocate (nKeep(0:lmax))
    if (Method_Mb_OrbBased_exactFieldFreeCut > 0.0_R64) then
      do l = 0, lmax
        nKeep(l) = count(blockEvals(:, l) < Method_Mb_OrbBased_exactFieldFreeCut)
      end do
      if (minval(nKeep) < 1) error stop "method.mb.orbBased.exactFieldFreeCut is below the lowest h0 eigenvalue"
      write (*, '(A,ES10.2,A,I0,A,I0,A,I0)') " spectral cut", Method_Mb_OrbBased_exactFieldFreeCut, &
        ": eigenvectors kept per l from ", minval(nKeep), " to ", maxval(nKeep), " of ", nRad
    else
      nKeep(:) = nRad
    end if

    !------------------------------------
    ! transformed (orbital, m) columns per l
    !------------------------------------

    ! keepQ(m + mCap + 1, iOrb), m running fastest
    nOS = Orbs_nOrbsInState
    allocate (blockColumns(0:lmax))
    do l = 0, lmax
      mCap = min(l, mmax)
      allocate (keepQ(2 * mCap + 1, nOS))
      keepQ(:, :) = .true.
      if (Method_Mb_OrbBased_definiteMQ) then
        keepQ(:, :) = .false.
        do iOrb = 1, nOS
          if (abs(Method_Mb_OrbBased_mOfOrb(iOrb)) <= mCap) &
            keepQ(Method_Mb_OrbBased_mOfOrb(iOrb) + mCap + 1, iOrb) = .true.
        end do
      end if
      blockColumns(l) % orb = pack(spread([(iOrb, iOrb=1, nOS)], 1, 2 * mCap + 1), keepQ)
      blockColumns(l) % m = pack(spread([(m, m=-mCap, mCap)], 2, nOS), keepQ)
      nColMax = max(nColMax, size(blockColumns(l) % orb))
      deallocate (keepQ)
    end do

    !------------------------------------
    ! optional GPU: eigenvectors stay on the device, staging buffers
    !------------------------------------

    gpuQ = Json_Get("method.mb.orbBased.exactFieldFreeGpuQ", .false.)
    if (gpuQ) then
      if (.not. CublasLib_AvailableQ()) &
        error stop "method.mb.orbBased.exactFieldFreeGpuQ: no CUDA device (or built without -DCODY_USE_CUDA=ON)"
      call CublasLib_Allocate(deviceEvecs, size(blockEvecs, kind=I64))
      call CublasLib_CopyToDevice(deviceEvecs, reshape(blockEvecs, [size(blockEvecs)]), 0_I64)
      blockSize = int(nRad, kind=I64) * 2 * nColMax
      allocate (hostBlocks(blockSize * (lmax + 1)))
      call CublasLib_Allocate(deviceX, blockSize * (lmax + 1))
      call CublasLib_Allocate(deviceY, blockSize * (lmax + 1))
      write (*, '(A,F8.1,A)') " exact field-free blocks on the GPU (", &
        real(size(blockEvecs, kind=I64) * 8, kind=R64) / 2.0_R64**20, " MiB)"
    end if

  end subroutine

end submodule
