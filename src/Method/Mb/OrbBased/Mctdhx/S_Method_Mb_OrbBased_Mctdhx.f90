! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased_Mctdhx) S_Method_Mb_OrbBased_Mctdhx

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Mctdhx_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Orbs
    use M_Coeffs
    use M_Method
    use M_Method_Mb_OrbBased

    call Say_Fabricate("method.manbody.orbbased.mctdhx")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Method_Setup => Setup
    Method_Mb_OrbBased_FillRdm1 => FillRdm1
    Method_Mb_OrbBased_FillRdm2 => FillRdm2
    Method_TimeDerivative => TimeDerivative
    Method_Mb_OrbBased_TimeDerivativeOrbsLin => TimeDerivativeOrbsLin
    Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin => TimeDerivativeCoeffsPlusOrbsNonLin
    Method_Mb_OrbBased_OrthonormalizeKeepingState => OrthonormalizeKeepingState
    if (Method_Mb_OrbBased_exactFieldFreeQ) then
      Method_PropagateExactPart => PropagateExactPart
      Method_ApplyExactPartPhi => ApplyExactPartPhi
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Orbs
    use M_OrbsInit
    use M_Coeffs
    use M_CoeffsInit
    use M_Method
    use M_Method_Mb_OrbBased

    integer(I32) :: nVals1, nVals2, nG, nOS, nC

    call Say_Setup("method.manbody.orbbased.mctdhx")
    call Method_Mb_OrbBased_Setup

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nC = Coeffs_nCoeffs

    nVals1 = nC
    nVals2 = nG * nOS

    allocate (Method_state(1:(nVals1 + nVals2)))

    call GetCoeffs(Coeffs_coeffs, Method_state)
    call CoeffsInit_Initialize(Coeffs_coeffs)

    call GetOrbs(Orbs_orbs, Method_state)
    call OrbsInit_Initialize(Orbs_orbs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm1(rdm1, state)
    use M_Coeffs

    complex(R64), intent(out), allocatable :: rdm1(:, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    complex(R64), contiguous, pointer :: coeffs(:)

    call GetCoeffs(coeffs, state)

    call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm2(rdm2, state)
    use M_Coeffs

    complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    complex(R64), contiguous, pointer :: coeffs(:)

    call GetCoeffs(coeffs, state)

    call Coeffs_ApplyH2FillRdm2(coeffs, rdm2_=rdm2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivative(dState, state, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    call FillTimeDerivative(dState, state, time, oneBodyQ=.true.)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Linear split-step part: the one-body operator on the orbitals only (the
  !> coefficients are constant), dOrbs = -i L orbs without projection,
  !> L = Method_Mb_OrbBased_ApplyLinOp (h, or only the gauge coupling if
  !> exactFieldFreeQ).
  subroutine TimeDerivativeOrbsLin(dState, state, time)
    use M_Method_Mb_OrbBased
    use M_Utils_Constants

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: dOrbs(:, :)

    dState(:) = 0.0_R64

    call GetOrbs(orbs, state)
    call GetOrbs(dOrbs, dState)

    call Method_Mb_OrbBased_ApplyLinOp(dOrbs, orbs, time)
    call Method_Mb_OrbBased_ApplySpectralCut(dOrbs)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Exact field-free propagation of the orbitals (the CI coefficients are
  !> constant in the linear split step).
  subroutine PropagateExactPart(state, t0, t1, imagTimeQ_)
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1
    logical, intent(in), optional :: imagTimeQ_

    complex(R64), contiguous, pointer :: orbs(:, :)
    logical :: imagTimeQ

    imagTimeQ = .false.
    if (present(imagTimeQ_)) imagTimeQ = imagTimeQ_

    call GetOrbs(orbs, state)
    call Method_Mb_OrbBased_PropagateFieldFree(orbs, t1 - t0, imagTimeQ)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Phi functions of the exact part: field-free h0 on the orbitals, zero
  !> generator (phi_m(0) = 1/m!) on the CI coefficients.
  subroutine ApplyExactPartPhi(res, vecs, dt)
    use M_Grid
    use M_Orbs
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous, target :: res(:)
    complex(R64), intent(in), contiguous, target :: vecs(:, :)
    real(R64), intent(in) :: dt

    complex(R64), contiguous, pointer :: coeffs(:), orbs(:, :)
    complex(R64), contiguous, pointer :: coeffsM(:), orbsM(:, :)
    complex(R64), allocatable :: orbsList(:, :, :)
    real(R64) :: invFactorial
    integer(I32) :: m

    call GetCoeffs(coeffs, res)
    call GetOrbs(orbs, res)
    allocate (orbsList(Grid_nPoints, Orbs_nOrbsInState, 0:size(vecs, 2) - 1))

    coeffs(:) = (0.0_R64, 0.0_R64)
    invFactorial = 1.0_R64
    do m = 0, size(vecs, 2) - 1
      if (m > 0) invFactorial = invFactorial / real(m, kind=R64)
      call GetCoeffs(coeffsM, vecs(:, m + 1))
      call GetOrbs(orbsM, vecs(:, m + 1))
      coeffs(:) = coeffs(:) + invFactorial * coeffsM(:)
      orbsList(:, :, m) = orbsM(:, :)
    end do

    call Method_Mb_OrbBased_ApplyFieldFreePhi(orbs, orbsList, dt)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivativeCoeffsPlusOrbsNonLin(dState, state, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    call FillTimeDerivative(dState, state, time, oneBodyQ=.false.)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Full time derivative (oneBodyQ), or the non-linear split-step part: only
  !> the interaction, since the linear part rotates the orbitals with the
  !> unprojected h at fixed coefficients (applying h1 here as well would count
  !> it twice); its orbital derivative is projected with the spectral cut.
  subroutine FillTimeDerivative(dState, state, time, oneBodyQ)
    use M_Coeffs
    use M_Orbs
    use M_Method_Mb_OrbBased
    use M_Utils_Constants

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time
    logical, intent(in)               :: oneBodyQ

    complex(R64), contiguous, pointer :: coeffs(:)
    complex(R64), contiguous, pointer :: dCoeffs(:)

    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: dOrbs(:, :)

    complex(R64), allocatable :: h1(:, :)
    complex(R64), allocatable :: h2(:, :, :, :)
    complex(R64), allocatable :: rdm1(:, :)
    complex(R64), allocatable :: invRdm1(:, :)
    complex(R64), allocatable :: rdm2(:, :, :, :)

    dState(:) = 0.0_R64

    call GetCoeffs(coeffs, state)
    call GetCoeffs(dCoeffs, dState)

    call GetOrbs(orbs, state)
    call GetOrbs(dOrbs, dState)

    call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1)
    call Coeffs_ApplyH2FillRdm2(coeffs, rdm2_=rdm2)

    if (oneBodyQ) call Method_Mb_OrbBased_ApplySingleBodyOp(dOrbs, orbs, time, h1_=h1)
    call Method_Mb_OrbBased_FillInverseRdm1(invRdm1, rdm1)
    call Method_Mb_OrbBased_ApplyCorrelationOp(dOrbs, orbs, invRdm1, rdm2, time, h2_=h2)

    if (oneBodyQ) call Coeffs_ApplyH1FillRdm1(coeffs, dCoeffs_=dCoeffs, h1_=h1)
    call Coeffs_ApplyH2FillRdm2(coeffs, dCoeffs_=dCoeffs, h2_=h2)

    call Coeffs_ProjectOnSubspace(dCoeffs, coeffs)
    call Orbs_ProjectOnSubspace(dOrbs, orbs)
    if (.not. oneBodyQ) call Method_Mb_OrbBased_ApplySpectralCut(dOrbs)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Loewdin orthonormalization phi <- phi A, A = S^(-1/2), of the stored
  !> orbitals (S = overlap matrix, zero between the blocks Orbs_blockOfOrb),
  !> with the CI transformation that keeps Psi = sum_I c_I Phi_I[phi] unchanged
  !> (exact in the full CI space of the active orbitals).
  !>
  !> With phi = (phi A) A^-1 the determinants transform with the one-body
  !> operator exp(K^) of K = -log(A) = log(S)/2, i.e. c <- exp(K^) c, where
  !> X^ = sum_kj X_kj a+_k a_j. The exponential is summed as a Taylor series of
  !> Coeffs_ApplyH1FillRdm1 applications (||log S|| is small).
  subroutine OrthonormalizeKeepingState(state)
    use M_Utils_LapackLib
    use M_Grid
    use M_Orbs
    use M_Coeffs
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous, target :: state(:)

    complex(R64), contiguous, pointer :: coeffs(:), orbs(:, :)
    real(R64), allocatable :: evals(:)
    complex(R64), allocatable :: evecs(:, :), s(:, :), kFull(:, :)
    complex(R64), allocatable :: term(:), termNew(:)
    integer(I32) :: nOS, i, j, k

    call GetCoeffs(coeffs, state)
    call GetOrbs(orbs, state)
    nOS = Orbs_nOrbsInState

    allocate (s(nOS, nOS))
    do j = 1, nOS
      do i = 1, nOS
        s(i, j) = 0.0_R64
        if (Orbs_blockOfOrb(i) .ne. Orbs_blockOfOrb(j)) cycle
        s(i, j) = Grid_InnerProduct(orbs(:, i), orbs(:, j))
      end do
    end do

    call LapackLib_DiagonalizeHermitianBlocks(evals, evecs, s, Orbs_blockOfOrb)

    ! Orbitals: phi <- phi S^(-1/2)
    orbs(:, :) = matmul(orbs, HermitianFunction(1.0_R64 / sqrt(evals)))

    ! K = log(S)/2 in the spin-orbital basis (restricted: same block for both spins)
    allocate (kFull(Method_Mb_OrbBased_nOrbsSum, Method_Mb_OrbBased_nOrbsSum))
    kFull(:, :) = 0.0_R64
    kFull(1:nOS, 1:nOS) = HermitianFunction(0.5_R64 * log(evals))
    if (Orbs_restrictedQ) kFull(nOS + 1:, nOS + 1:) = kFull(1:nOS, 1:nOS)

    ! c <- exp(K^) c
    term = coeffs
    allocate (termNew, mold=coeffs)
    do k = 1, 60
      termNew(:) = 0.0_R64
      call Coeffs_ApplyH1FillRdm1(term, dCoeffs_=termNew, h1_=kFull)
      term = termNew / k
      coeffs = coeffs + term
      if (maxval(abs(term)) < 1.0e-16_R64) exit
    end do

  contains

    !> U diag(f) U^H with the eigenvectors U of S.
    function HermitianFunction(f) result(res)
      real(R64), intent(in), contiguous :: f(:)
      complex(R64), allocatable :: res(:, :)

      res = matmul(evecs, spread(f, 2, size(f)) * conjg(transpose(evecs)))

    end function

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GetCoeffs(coeffs, state)
    use M_Coeffs
    complex(R64), pointer, intent(out) :: coeffs(:)
    complex(R64), intent(in), contiguous, target :: state(:)

    coeffs(1:Coeffs_nCoeffs) => state(1:Coeffs_nCoeffs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GetOrbs(orbs, state)
    use M_Coeffs
    use M_Grid
    use M_Orbs
    use M_Method_Mb_OrbBased
    complex(R64), pointer, intent(out) :: orbs(:, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    orbs(1:Grid_nPoints, 1:Orbs_nOrbsInState) => state(Coeffs_nCoeffs + 1:)

  end subroutine

end submodule
