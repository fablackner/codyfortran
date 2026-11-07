! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief DIIS (Pulay/Anderson) mixing implementation submodule.
!>
!> @details
!> See M_Mixing_Diis. The history stores, per pair, the residual
!> r_i = xRaw_i − xMixed_i and the damped candidate y_i = xMixed_i + α·r_i,
!> so that the step Σ c_i·(xMixed_i + α·r_i) = Y·c is one matrix-vector
!> product. The real Pulay matrix B_ij = Re⟨r_i|r_j⟩ is kept across calls:
!> a new pair adds one row (one matrix-vector product with the stored
!> residuals), so a step costs O(m·dim) instead of O(m²·dim).
submodule(M_Mixing_Diis) S_Mixing_Diis

  implicit none

  !> Relative eigenvalue cutoff of the pseudo-inverse of the bordered Pulay
  !> matrix.
  real(R64), parameter :: eigCutoff = 1e-12_R64

  !> Restart threshold on Σ|c_i|: beyond this the extrapolation is considered
  !> untrustworthy (near-dependent residuals) and the history is dropped.
  real(R64), parameter :: coeffLimit = 1e4_R64

  !> Residual growth factor that triggers a history restart.
  real(R64), parameter :: growthLimit = 3.0_R64

  !=============================================================================
  ! configuration (read from JSON in Fabricate)
  !=============================================================================

  !> Maximum number of stored (candidate, residual) pairs.
  integer(I32) :: nHistory = 8

  !> Start gate: plain damped steps while the residual norm is above this
  !> factor times the first residual norm (0: no gate).
  real(R64) :: startThreshold = 0.1_R64

  !> Damping parameter (0 < α ≤ 1).
  real(R64) :: alpha = 1.0_R64

  !> 0: silent; 1: one line per step (residual norm, pairs, restarts).
  integer(I32) :: printLevel = 0

  !=============================================================================
  ! mixer state (module-level, lazily allocated)
  !=============================================================================

  !> Length of the mixed vectors (0: not allocated).
  integer(I32) :: dim = 0

  !> Number of pairs currently stored (slots 1..nStored).
  integer(I32) :: nStored = 0

  !> Slot of the most recently stored pair.
  integer(I32) :: iNewest = 0

  !> Steps since the last reset.
  integer(I32) :: nSteps = 0

  !> Residual norm of the first step since the last reset (start-gate reference).
  real(R64) :: rNormFirst = -1.0_R64

  !> Residual norm of the previous step (growth detection).
  real(R64) :: rNormPrev = huge(1.0_R64)

  !> Damped candidates y_i = xMixed_i + α·r_i, dimension(dim, nHistory).
  complex(R64), allocatable :: yHistory(:, :)

  !> Residuals r_i = xRaw_i − xMixed_i, dimension(dim, nHistory).
  complex(R64), allocatable :: rHistory(:, :)

  !> Pulay matrix B_ij = Re⟨r_i|r_j⟩ of the stored slots, dimension(nHistory, nHistory).
  real(R64), allocatable :: pulayMat(:, :)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Mixing_Diis_Fabricate()
    use M_Utils_Say
    use M_Utils_Json
    use M_Mixing

    call Say_Fabricate("mixing.diis")

    !------------------------------------
    ! read configuration parameters
    !------------------------------------

    nHistory = Json_Get("nHistory", 8, path_="mixing.diis")
    if (nHistory < 1) error stop "mixing.diis.nHistory must be >= 1"

    startThreshold = Json_Get("startThreshold", 0.1_R64, path_="mixing.diis")
    if (startThreshold < 0.0_R64) error stop "mixing.diis.startThreshold must be >= 0 (0: no start gate)"

    alpha = Json_Get("alpha", 1.0_R64, path_="mixing.diis")
    if (.not. (alpha > 0.0_R64 .and. alpha <= 1.0_R64)) error stop "mixing.diis.alpha must be in (0, 1]"

    printLevel = Json_Get("printLevel", 0, path_="mixing.diis")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! a new mixer starts empty: buffers are (re)allocated at the first Mix
    dim = 0
    call ResetMixer

    Mixing_Mix => Mix
    Mixing_Reset => ResetMixer

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Drop the history and the start-gate reference.
  subroutine ResetMixer

    nStored = 0
    iNewest = 0
    nSteps = 0
    rNormFirst = -1.0_R64
    rNormPrev = huge(1.0_R64)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Allocate the history for vectors of length newDim (a new length
  !> resets the mixer).
  subroutine EnsureAllocated(newDim)
    integer(I32), intent(in) :: newDim

    if (dim .eq. newDim) return

    dim = newDim
    call ResetMixer

    if (allocated(yHistory)) deallocate (yHistory)
    if (allocated(rHistory)) deallocate (rHistory)
    if (allocated(pulayMat)) deallocate (pulayMat)
    allocate (yHistory(dim, nHistory), rHistory(dim, nHistory), pulayMat(nHistory, nHistory))
    pulayMat(:, :) = 0.0_R64

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief One DIIS step: store the pair of (xMixed, xRaw), extrapolate,
  !> damped update (see M_Mixing_Diis).
  subroutine Mix(xMixed, xRaw)
    use M_Utils_BlasLib

    complex(R64), intent(inout), contiguous :: xMixed(:)
    complex(R64), intent(in), contiguous :: xRaw(:)

    complex(R64), allocatable :: overlaps(:)
    real(R64), allocatable :: coeffs(:)
    real(R64) :: rNorm
    character(len=:), allocatable :: note
    logical :: okQ
    integer(I32) :: m

    call EnsureAllocated(size(xMixed))
    nSteps = nSteps + 1

    ! Store the pair in the circular buffer: residual and damped candidate
    iNewest = mod(iNewest, nHistory) + 1
    nStored = min(nStored + 1, nHistory)
    rHistory(:, iNewest) = xRaw - xMixed
    yHistory(:, iNewest) = xMixed + alpha * rHistory(:, iNewest)

    ! New row of the Pulay matrix: ⟨r_i|r_newest⟩ for all stored slots (the
    ! slots 1..nStored are filled)
    allocate (overlaps(nStored))
    call BlasLib_CalcMatOnVec(overlaps, rHistory(:, 1:nStored), rHistory(:, iNewest), 'C')
    pulayMat(1:nStored, iNewest) = real(overlaps, kind=R64)
    pulayMat(iNewest, 1:nStored) = real(overlaps, kind=R64)
    rNorm = sqrt(max(pulayMat(iNewest, iNewest), 0.0_R64))

    ! Start gate: extrapolation trusts a linear model of the fixed-point map,
    ! which holds only near convergence: plain damped steps (newest pair only)
    ! while the residual is above startThreshold × the first residual
    ! (startThreshold 0: no gate)
    note = ""
    if (rNormFirst < 0.0_R64) rNormFirst = rNorm
    if (0.0_R64 < startThreshold .and. startThreshold * rNormFirst < rNorm) then
      call RestartWithNewestPair
    else if (growthLimit * rNormPrev < rNorm) then
      ! the residual jumped: the history describes another region of the map
      call RestartWithNewestPair
      note = " (restart: residual growth)"
    end if
    rNormPrev = rNorm

    ! Pulay coefficients: minimize ‖Σ c_i r_i‖² subject to Σ c_i = 1
    m = nStored
    allocate (coeffs(m))
    coeffs(:) = 1.0_R64
    if (1 < m) then
      call PulayCoefficients(coeffs, okQ, pulayMat(1:m, 1:m))
      ! Safeguard: no solution or large coefficients mean (nearly) linearly
      ! dependent residuals and an untrustworthy linear model
      if (.not. okQ) then
        call RestartWithNewestPair
        m = 1
        coeffs = [1.0_R64]
        note = " (restart: dependent residuals)"
      end if
    end if

    if (0 < printLevel) &
      write (*, '(1X, A, I4, A, ES10.3, A, I2, A)') "mixing.diis: step ", nSteps, ", |r| ", rNorm, ", pairs ", m, note

    ! Damped step along the extrapolated pair: xMixed = Σ c_i y_i (with one
    ! stored pair: (1 − α)·xMixed + α·xRaw)
    call BlasLib_CalcMatOnVec(xMixed, yHistory(:, 1:m), cmplx(coeffs, 0.0_R64, kind=R64), 'N')

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Real Pulay coefficients of the stored pairs.
  !>
  !> @details
  !> Minimizes cᵀBc subject to Σc = 1 (B = Re⟨r_i|r_j⟩, scaled to unit
  !> maximum diagonal: changes only the Lagrange multiplier) through the
  !> bordered system
  !>
  !>   [ B   1 ] [c]   [0]
  !>   [ 1ᵀ  0 ] [λ] = [1]
  !>
  !> solved with the pseudo-inverse of its symmetric eigen-decomposition
  !> (eigenvalues below eigCutoff × the largest |eigenvalue| dropped: the
  !> SVD pseudo-inverse of a symmetric matrix), which stays well-behaved for
  !> (nearly) linearly dependent residuals close to convergence. Real
  !> coefficients keep the mixed iterate in the real-affine hull of the
  !> history: complex ones would rotate essentially real quantities
  !> (densities, potentials) into their imaginary parts. okQ is false if the
  !> coefficient sum degenerates (|Σc| < 1e-8 before renormalization) or
  !> Σ|c| exceeds coeffLimit.
  subroutine PulayCoefficients(coeffs, okQ, b)
    use M_Utils_LapackLib

    real(R64), intent(out), contiguous :: coeffs(:)
    logical, intent(out) :: okQ
    real(R64), intent(in), contiguous :: b(:, :)

    real(R64), allocatable :: bordered(:, :), evals(:), evecs(:, :), z(:)
    real(R64) :: maxDiag, coeffSum
    integer(I32) :: k, m

    m = size(b, 1)
    maxDiag = 0.0_R64
    do k = 1, m
      maxDiag = max(maxDiag, b(k, k))
    end do

    allocate (bordered(m + 1, m + 1))
    bordered(1:m, 1:m) = b
    if (0.0_R64 < maxDiag) bordered(1:m, 1:m) = b / maxDiag
    bordered(1:m, m + 1) = 1.0_R64
    bordered(m + 1, 1:m) = 1.0_R64
    bordered(m + 1, m + 1) = 0.0_R64

    ! [c; λ] = V Λ⁺ Vᵀ e_{m+1}
    call LapackLib_DiagonalizeGeneric(evals, evecs, bordered, .true.)
    z = evecs(m + 1, :)
    do k = 1, m + 1
      if (eigCutoff * maxval(abs(evals)) < abs(evals(k))) then
        z(k) = z(k) / evals(k)
      else
        z(k) = 0.0_R64
      end if
    end do
    coeffs(:) = matmul(evecs(1:m, :), z)

    ! renormalize the constraint (exact solves satisfy it; the cutoff can
    ! introduce small deviations)
    coeffSum = sum(coeffs)
    okQ = 1e-8_R64 < abs(coeffSum)
    if (okQ) then
      coeffs(:) = coeffs(:) / coeffSum
      okQ = sum(abs(coeffs)) < coeffLimit
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Drop the history except for the most recent pair (moved to slot 1).
  subroutine RestartWithNewestPair

    if (1 < iNewest) then
      yHistory(:, 1) = yHistory(:, iNewest)
      rHistory(:, 1) = rHistory(:, iNewest)
      pulayMat(1, 1) = pulayMat(iNewest, iNewest)
    end if
    nStored = 1
    iNewest = 1

  end subroutine

end submodule
