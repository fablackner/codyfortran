! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief DIIS (Pulay/Anderson) mixing implementation submodule.
submodule(M_Mixing_Diis) S_Mixing_Diis

  implicit none

  !> Relative singular-value cutoff for the pseudo-inverse of the Pulay matrix.
  real(R64), parameter :: svdCutoff = 1e-12_R64

  !> Restart threshold on Σ|c_i|: beyond this the extrapolation is considered
  !> untrustworthy (near-dependent residuals) and the history is dropped.
  real(R64), parameter :: coeffLimit = 1e4_R64

  !> Residual growth factor that triggers a history restart.
  real(R64), parameter :: growthLimit = 3.0_R64

  !=============================================================================
  ! configuration (read from JSON in Fabricate)
  !=============================================================================

  !> Maximum number of stored (iterate, residual) pairs.
  integer(I32) :: nHistory = 8

  !> Extrapolation starts once the residual norm has dropped below this factor
  !> times the first residual norm; until then plain damped steps are taken.
  real(R64) :: startThreshold = 0.1_R64

  !> Damping parameter (0 < α ≤ 1). Read from JSON during fabrication.
  real(R64) :: alpha = 1.0_R64

  !=============================================================================
  ! mixer state (module-level, lazily allocated)
  !=============================================================================

  !> Length of the mixed vectors (set on first call / dimension change).
  integer(I32) :: dim = 0

  !> Number of (iterate, residual) pairs currently stored.
  integer(I32) :: nStored = 0

  !> Circular-buffer slot of the most recently stored pair.
  integer(I32) :: iNewest = 0

  !> History of accepted iterates x_i, dimension(dim, nHistory).
  complex(R64), allocatable :: xHistory(:, :)

  !> History of residuals r_i = xRaw_i − xMixed_i, dimension(dim, nHistory).
  complex(R64), allocatable :: rHistory(:, :)

  !> Residual norm of the first step since setup/reset (DIIS gating reference).
  real(R64) :: rNormFirst = -1.0_R64

  !> Residual norm of the previous step (DIIS growth detection).
  real(R64) :: rNormPrev = huge(1.0_R64)

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

    alpha = Json_Get("alpha", 1.0_R64, path_="mixing.diis")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Mixing_Mix => Mix

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Ensure the history buffers are allocated for the current vector dimension.
  !>
  !> @details
  !> Called at the start of each Mix step. If the dimension has changed
  !> (or this is the first call), the history is (re)allocated and reset.
  subroutine EnsureAllocated(newDim)
    integer(I32), intent(in) :: newDim

    if (dim .eq. newDim) return

    dim = newDim
    nStored = 0
    iNewest = 0
    rNormFirst = -1.0_R64
    rNormPrev = huge(1.0_R64)

    if (allocated(xHistory)) deallocate (xHistory)
    if (allocated(rHistory)) deallocate (rHistory)
    allocate (xHistory(dim, nHistory))
    allocate (rHistory(dim, nHistory))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief One DIIS step: store (xMixed, xRaw−xMixed), extrapolate, damped update.
  !>
  !> @details
  !> The Pulay coefficients are obtained from the bordered system solved with
  !> an SVD pseudo-inverse; singular values below `svdCutoff` relative to the
  !> largest are dropped, which keeps the extrapolation stable when residuals
  !> become (nearly) linearly dependent close to convergence. If the
  !> coefficient sum degenerates, the step falls back to plain linear mixing.
  subroutine Mix(xMixed, xRaw)
    use M_Utils_LapackLib

    complex(R64), intent(inout), contiguous :: xMixed(:)
    complex(R64), intent(in), contiguous :: xRaw(:)

    complex(R64), allocatable :: pulayMat(:, :), rhs(:), coeffs(:), z(:)
    complex(R64), allocatable :: u(:, :), vt(:, :)
    real(R64), allocatable :: s(:)
    complex(R64) :: coeffSum
    real(R64) :: maxDiag, rNormNew
    integer(I32) :: m, i, j, k

    call EnsureAllocated(size(xMixed))

    ! Store the pair (xMixed, r = xRaw − xMixed) in the circular history buffer
    iNewest = mod(iNewest, nHistory) + 1
    xHistory(:, iNewest) = xMixed
    rHistory(:, iNewest) = xRaw - xMixed
    nStored = min(nStored + 1, nHistory)

    ! Gating: extrapolation trusts a linear model of the fixed-point map,
    ! which only holds near convergence. Until the residual norm has dropped
    ! below startThreshold × (first residual norm), and whenever it grows
    ! sharply, keep only the newest pair — the solve below then degenerates to
    ! a plain damped linear step.
    rNormNew = sqrt(real(dot_product(rHistory(:, iNewest), &
                                     rHistory(:, iNewest)), R64))
    if (rNormFirst < 0.0_R64) rNormFirst = rNormNew

    if (startThreshold * rNormFirst < rNormNew) then
      call RestartWithNewestPair
    else if (growthLimit * rNormPrev < rNormNew) then
      call RestartWithNewestPair
    end if
    rNormPrev = rNormNew

    m = nStored

    ! Bordered Pulay system: minimize ‖Σ c_i·r_i‖² subject to Σ c_i = 1.
    ! The coefficients are restricted to be real (only Re⟨r_i|r_j⟩ enters), so
    ! the mixed iterate stays in the real-affine hull of the history — complex
    ! coefficients would rotate essentially-real quantities (densities,
    ! potentials) into their imaginary parts.
    allocate (pulayMat(m + 1, m + 1))
    allocate (rhs(m + 1))

    do j = 1, m
      do i = 1, m
        pulayMat(i, j) = real(dot_product(rHistory(:, i), rHistory(:, j)), R64)
      end do
    end do

    ! Scale the residual block to O(1): changes only the Lagrange multiplier,
    ! not the coefficients, but keeps the SVD cutoff meaningful when the
    ! residual norms are far from unity
    maxDiag = 0.0_R64
    do i = 1, m
      maxDiag = max(maxDiag, real(pulayMat(i, i), R64))
    end do
    if (0.0_R64 < maxDiag) pulayMat(1:m, 1:m) = pulayMat(1:m, 1:m) / maxDiag

    pulayMat(1:m, m + 1) = (1.0_R64, 0.0_R64)
    pulayMat(m + 1, 1:m) = (1.0_R64, 0.0_R64)
    pulayMat(m + 1, m + 1) = (0.0_R64, 0.0_R64)

    rhs = (0.0_R64, 0.0_R64)
    rhs(m + 1) = (1.0_R64, 0.0_R64)

    ! Solve via SVD pseudo-inverse: coeffs = V·Σ⁺·U†·rhs
    call LapackLib_Svd(u, s, vt, pulayMat)

    allocate (z(m + 1))
    z = matmul(conjg(transpose(u)), rhs)
    do k = 1, m + 1
      if (svdCutoff * s(1) < s(k)) then
        z(k) = z(k) / s(k)
      else
        z(k) = (0.0_R64, 0.0_R64)
      end if
    end do
    coeffs = matmul(conjg(transpose(vt)), z)

    ! Renormalize the constraint Σ c_i = 1 (exact solves already satisfy it;
    ! the pseudo-inverse cutoff can introduce small deviations)
    coeffSum = sum(coeffs(1:m))
    if (abs(coeffSum) < 1e-8_R64) then
      call RestartWithNewestPair
      xMixed = (1.0_R64 - alpha) * xMixed + alpha * xRaw
      return
    end if
    coeffs(1:m) = coeffs(1:m) / coeffSum

    ! Safeguard: large coefficients mean (nearly) linearly dependent residuals
    ! and an untrustworthy linear model — restart the history and take a plain
    ! damped step instead of extrapolating wildly
    if (coeffLimit < sum(abs(coeffs(1:m)))) then
      call RestartWithNewestPair
      xMixed = (1.0_R64 - alpha) * xMixed + alpha * xRaw
      return
    end if

    ! Damped step along the extrapolated iterate/residual pair
    xMixed = (0.0_R64, 0.0_R64)
    do i = 1, m
      xMixed = xMixed + coeffs(i) * (xHistory(:, i) + alpha * rHistory(:, i))
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Drop the history except for the most recent (iterate, residual) pair.
  subroutine RestartWithNewestPair

    if (1 < iNewest) then
      xHistory(:, 1) = xHistory(:, iNewest)
      rHistory(:, 1) = rHistory(:, iNewest)
    end if
    nStored = 1
    iNewest = 1

  end subroutine

end submodule
