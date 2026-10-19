! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased) S_Method_Mb_OrbBased_FieldFree
  !-----------------------------------------------------------------------------
  ! Exact propagation with the field-free one-body Hamiltonian h0 = T + V_ext
  ! on Ylm grids (method.mb.orbBased.exactFieldFreeQ).
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

  !> GPU variant (method.mb.orbBased.exactFieldFreeGpuQ): eigenvectors on the
  !> device, device work arrays and the host staging buffer.
  logical :: gpuQ = .false.
  type(c_ptr) :: deviceEvecs = c_null_ptr
  type(c_ptr) :: deviceX = c_null_ptr
  type(c_ptr) :: deviceY = c_null_ptr
  real(R64), allocatable, target :: hostBlocks(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_PropagateFieldFree(orbs, dt, imagTimeQ)
    use M_Utils_BlasLib, only: DGEMM
    use M_Grid_Ylm

    complex(R64), intent(inout), contiguous :: orbs(:, :)
    real(R64), intent(in) :: dt
    logical, intent(in) :: imagTimeQ

    integer(I32) :: nRad, l, nCol
    real(R64), allocatable :: x(:, :), y(:, :)

    if (.not. allocated(blockEvecs)) call BuildBlocks

    if (gpuQ) then
      call PropagateFieldFreeGpu(orbs, dt, imagTimeQ)
      return
    end if

    nRad = Grid_Ylm_nRadial

    !$omp parallel do default(shared) schedule(dynamic) private(l, nCol, x, y)
    do l = 0, Grid_Ylm_lmax
      nCol = NColumns(l, size(orbs, 2))
      allocate (x(nRad, 2 * nCol), y(nRad, 2 * nCol))

      ! Eigenbasis coefficients U^T x, phases, back transformation U y
      call GatherBlock(x, orbs, l)
      call DGEMM('T', 'N', nRad, 2 * nCol, nRad, 1.0_R64, blockEvecs(:, :, l), nRad, x, nRad, 0.0_R64, y, nRad)
      call ApplyPhases(y, l, dt, imagTimeQ)
      call DGEMM('N', 'N', nRad, 2 * nCol, nRad, 1.0_R64, blockEvecs(:, :, l), nRad, y, nRad, 0.0_R64, x, nRad)
      call ScatterBlock(orbs, x, l)

      deallocate (x, y)
    end do
    !$omp end parallel do

  end subroutine

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

    integer(I32) :: nRad, lmax, l, nCol, nColMax
    integer(I64) :: blockSize, offsetU, offsetX
    real(R64), contiguous, pointer :: blocks(:, :, :)

    nRad = Grid_Ylm_nRadial
    lmax = Grid_Ylm_lmax
    nColMax = NColumns(lmax, size(orbs, 2))
    blockSize = int(nRad, kind=I64) * 2 * nColMax

    ! Host staging buffer and device work arrays (2 nColMax columns per l)
    if (.not. allocated(hostBlocks)) then
      allocate (hostBlocks(blockSize * (lmax + 1)))
      call CublasLib_Allocate(deviceX, blockSize * (lmax + 1))
      call CublasLib_Allocate(deviceY, blockSize * (lmax + 1))
    end if
    if (size(hostBlocks, kind=I64) .ne. blockSize * (lmax + 1)) &
      error stop "Method_Mb_OrbBased_PropagateFieldFree: orbital count changed between GPU calls"

    blocks(1:nRad, 1:2 * nColMax, 0:lmax) => hostBlocks

    !$omp parallel do default(shared) schedule(dynamic) private(l, nCol)
    do l = 0, lmax
      nCol = NColumns(l, size(orbs, 2))
      call GatherBlock(blocks(:, 1:2 * nCol, l), orbs, l)
    end do
    !$omp end parallel do

    call CublasLib_CopyToDevice(deviceX, hostBlocks, 0_I64)
    do l = 0, lmax
      nCol = NColumns(l, size(orbs, 2))
      offsetU = int(nRad, kind=I64) * nRad * l
      offsetX = blockSize * l
      call CublasLib_Dgemm('T', 'N', nRad, 2 * nCol, nRad, 1.0_R64, deviceEvecs, offsetU, nRad, &
                           deviceX, offsetX, nRad, 0.0_R64, deviceY, offsetX, nRad)
    end do
    call CublasLib_CopyToHost(hostBlocks, deviceY, 0_I64)

    !$omp parallel do default(shared) schedule(dynamic) private(l, nCol)
    do l = 0, lmax
      nCol = NColumns(l, size(orbs, 2))
      call ApplyPhases(blocks(:, 1:2 * nCol, l), l, dt, imagTimeQ)
    end do
    !$omp end parallel do

    call CublasLib_CopyToDevice(deviceY, hostBlocks, 0_I64)
    do l = 0, lmax
      nCol = NColumns(l, size(orbs, 2))
      offsetU = int(nRad, kind=I64) * nRad * l
      offsetX = blockSize * l
      call CublasLib_Dgemm('N', 'N', nRad, 2 * nCol, nRad, 1.0_R64, deviceEvecs, offsetU, nRad, &
                           deviceY, offsetX, nRad, 0.0_R64, deviceX, offsetX, nRad)
    end do
    call CublasLib_CopyToHost(hostBlocks, deviceX, 0_I64)

    !$omp parallel do default(shared) schedule(dynamic) private(l, nCol)
    do l = 0, lmax
      nCol = NColumns(l, size(orbs, 2))
      call ScatterBlock(orbs, blocks(:, 1:2 * nCol, l), l)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Number of (orbital, m) columns of the block of l.
  pure function NColumns(l, nOS) result(res)
    use M_Grid_Ylm

    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: nOS
    integer(I32) :: res

    res = (2 * min(l, Grid_Ylm_mmax) + 1) * nOS

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Gathers w^{1/2} f_lm of all orbitals and m of the block of l: real parts
  !> in columns 1:nCol, imaginary parts in nCol+1:2 nCol (one GEMM per
  !> transformation).
  subroutine GatherBlock(x, orbs, l)
    use M_Grid_Ylm

    real(R64), intent(out), contiguous :: x(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    integer(I32), intent(in) :: l

    integer(I32) :: iOrb, m, mCap, iCol, nCol
    complex(R64), allocatable :: fLm(:)

    allocate (fLm(Grid_Ylm_nRadial))
    mCap = min(l, Grid_Ylm_mmax)
    nCol = size(x, 2) / 2

    iCol = 0
    do iOrb = 1, size(orbs, 2)
      do m = -mCap, mCap
        iCol = iCol + 1
        call Grid_Ylm_GetLmComponent(fLm, l, m, orbs(:, iOrb))
        x(:, iCol) = sqrtWeights(:) * real(fLm(:), kind=R64)
        x(:, nCol + iCol) = sqrtWeights(:) * aimag(fLm(:))
      end do
    end do

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
    if (imagTimeQ) then
      phase(:) = cmplx(exp(-blockEvals(:, l) * dt), 0.0_R64, kind=R64)
    else
      phase(:) = cmplx(cos(blockEvals(:, l) * dt), -sin(blockEvals(:, l) * dt), kind=R64)
    end if

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
  !> Scatters w^{-1/2} x of the block of l back into the channels.
  subroutine ScatterBlock(orbs, x, l)
    use M_Grid_Ylm

    complex(R64), intent(inout), contiguous :: orbs(:, :)
    real(R64), intent(in), contiguous :: x(:, :)
    integer(I32), intent(in) :: l

    integer(I32) :: iOrb, m, mCap, iCol, nCol
    complex(R64), allocatable :: fLm(:)

    allocate (fLm(Grid_Ylm_nRadial))
    mCap = min(l, Grid_Ylm_mmax)
    nCol = size(x, 2) / 2

    iCol = 0
    do iOrb = 1, size(orbs, 2)
      do m = -mCap, mCap
        iCol = iCol + 1
        fLm(:) = cmplx(x(:, iCol), x(:, nCol + iCol), kind=R64) / sqrtWeights(:)
        call Grid_Ylm_SetLmComponent(orbs(:, iOrb), l, m, fLm)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine BuildBlocks
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

    real(R64), parameter :: tolerance = 1.0e-10_R64

    integer(I32) :: nRad, nCh, nG, lmax, mmax, iCh, i, l, m
    logical :: potentialQ
    real(R64) :: hScale, devM, devImag, devSym
    complex(R64), allocatable :: probe(:), hProbe(:), tmp(:), hLm(:), externalPotential(:)
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
    nCh = Grid_Ylm_NChannels(lmax, mmax)

    allocate (sqrtWeights(nRad))
    sqrtWeights(:) = sqrt(Grid_Ylm_radialWeights(:))

    allocate (probe(nG), hProbe(nG), tmp(nG), hLm(nRad))
    allocate (blockEvecs(nRad, nRad, 0:lmax), blockEvals(nRad, 0:lmax))

    if (potentialQ) call SysPotential_FillExternalPotential(externalPotential, 0.0_R64)

    devM = 0.0_R64
    devImag = 0.0_R64

    do i = 1, nRad
      probe(:) = (0.0_R64, 0.0_R64)
      do iCh = 1, nCh
        probe((iCh - 1) * nRad + i) = (1.0_R64, 0.0_R64)
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

    !$omp parallel do default(shared) schedule(dynamic) private(l, hs, evals, evecs)
    do l = 0, lmax
      hs = 0.5_R64 * (blockEvecs(:, :, l) + transpose(blockEvecs(:, :, l)))
      call LapackLib_DiagonalizeGeneric(evals, evecs, hs, .true.)
      blockEvals(:, l) = evals(:)
      blockEvecs(:, :, l) = evecs(:, :)
    end do
    !$omp end parallel do

    write (*, '(A,I0,A,ES12.4,A,ES12.4)') " exact field-free blocks diagonalized for l <= ", lmax, &
      ": eigenvalues from", minval(blockEvals), " to", maxval(blockEvals)

    !------------------------------------
    ! optional GPU: eigenvectors stay on the device
    !------------------------------------

    gpuQ = Json_Get("method.mb.orbBased.exactFieldFreeGpuQ", .false.)
    if (gpuQ) then
      if (.not. CublasLib_AvailableQ()) &
        error stop "method.mb.orbBased.exactFieldFreeGpuQ: no CUDA device (or built without -DCODY_USE_CUDA=ON)"
      call CublasLib_Allocate(deviceEvecs, size(blockEvecs, kind=I64))
      call CublasLib_CopyToDevice(deviceEvecs, reshape(blockEvecs, [size(blockEvecs)]), 0_I64)
      write (*, '(A,F8.1,A)') " exact field-free blocks on the GPU (", &
        real(size(blockEvecs, kind=I64) * 8, kind=R64) / 2.0_R64**20, " MiB)"
    end if

  end subroutine

end submodule
