!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Helpers shared by the beryllium programs (included after `contains`).
!
! - EmbedOrbitals:        load orbXX_YY.in written on a grid with a smaller
!                         lmax (same radial grid, mmax = 0)
! - TransformOrbitalsKeepingState:
!                         in-space orbital transformation phi -> phi A with the
!                         compensating CI transformation, so that the
!                         many-body state is unchanged
! - LowdinKeepingState:   symmetric orthonormalization of the orbitals with the
!                         compensating CI transformation
! - MetricWeights:        FEDVR metric of a packed state (SIL integrator)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fill the orbitals from orbXX_YY.in files written on a Ylm grid with the same
!> radial grid but a smaller lmax (lmaxIn). With mmax = 0 the channel layout
!> is l-major, so the stored orbital is a prefix of the current one and the
!> missing high-l channels are zero.
subroutine EmbedOrbitals(orbs, lmaxIn)
  use M_Utils_DataStorage
  use M_Grid_Ylm
  use M_Method_Mb
  use M_Method_Mb_OrbBased
  use M_Orbs

  complex(R64), intent(inout), contiguous :: orbs(:, :)
  integer(I32), intent(in) :: lmaxIn

  complex(R64), allocatable :: orbIn(:)
  character(len=256) :: filename
  integer(I32) :: nBt, ibt, index, i1, nIn

  if (Grid_Ylm_mmax .ne. 0) error stop "EmbedOrbitals requires grid.ylm.mmax = 0"
  if (lmaxIn > Grid_Ylm_lmax) error stop "EmbedOrbitals: lmaxIn exceeds grid.ylm.lmax"

  nIn = (lmaxIn + 1) * Grid_Ylm_nRadial
  allocate (orbIn(nIn))

  nBt = Method_Mb_nBodyTypes
  if (Orbs_restrictedQ) nBt = 1

  orbs(:, :) = 0.0_R64

  i1 = 0
  do ibt = 1, nBt
    do index = 1, Method_Mb_OrbBased_nOrbs(ibt)
      i1 = i1 + 1
      write (filename, '("orb",I2.2,"_",I2.2,".in")') ibt, index
      call LoadData(orbIn, trim(filename), storage_size(orbIn))
      orbs(1:nIn, i1) = orbIn(:)
    end do
  end do

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Replace the orbitals by phi A for a Hermitian positive-definite in-space
!> matrix A (one block per stored orbital set) and transform the CI vector so
!> that Psi = sum_I c_I Phi_I[phi] is unchanged.
!>
!> With phi = (phi A) A^-1 the determinants transform with the one-body
!> operator exp(K^) of K = log(A^-1) = -log(A), i.e. c <- exp(-log(A)^) c,
!> where X^ = sum_kj X_kj a+_k a_j. The exponential is summed as a Taylor
!> series of Coeffs_ApplyH1FillRdm1 applications (||log A|| is small).
subroutine TransformOrbitalsKeepingState(coeffs, orbs, a)
  use M_Utils_LapackLib
  use M_Grid
  use M_Coeffs
  use M_Orbs
  use M_Method_Mb_OrbBased

  complex(R64), intent(inout), contiguous :: coeffs(:)
  complex(R64), intent(inout), contiguous :: orbs(:, :)
  complex(R64), intent(in), contiguous    :: a(:, :)

  real(R64), allocatable :: evals(:)
  complex(R64), allocatable :: evecs(:, :), logA(:, :), kFull(:, :)
  complex(R64), allocatable :: term(:), termNew(:)
  integer(I32) :: nOS, nO, k

  nOS = Orbs_nOrbsInState
  nO = Method_Mb_OrbBased_nOrbsSum

  ! Orbitals: phi <- phi A
  orbs(:, 1:nOS) = matmul(orbs(:, 1:nOS), a)

  ! K = -log(A) in the spin-orbital basis (restricted: same block for both spins)
  call LapackLib_DiagonalizeGeneric(evals, evecs, a, .true.)
  if (minval(evals) <= 0.0_R64) error stop "TransformOrbitalsKeepingState: A is not positive definite"
  logA = matmul(evecs, spread(log(evals), 2, size(evals)) * conjg(transpose(evecs)))

  allocate (kFull(nO, nO))
  kFull(:, :) = 0.0_R64
  kFull(1:nOS, 1:nOS) = -logA
  if (Orbs_restrictedQ) kFull(nOS + 1:, nOS + 1:) = -logA

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

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Loewdin orthonormalization phi <- phi S^(-1/2) of the stored orbitals
!> (S = overlap matrix, cross body-type blocks are zero), keeping the many-body
!> state unchanged through the compensating CI transformation.
subroutine LowdinKeepingState(coeffs, orbs)
  use M_Utils_LapackLib
  use M_Grid
  use M_Orbs
  use M_Method_Mb_OrbBased

  complex(R64), intent(inout), contiguous :: coeffs(:)
  complex(R64), intent(inout), contiguous :: orbs(:, :)

  real(R64), allocatable :: evals(:)
  complex(R64), allocatable :: evecs(:, :), s(:, :), sInvSqrt(:, :)
  integer(I32) :: nOS, i, j

  nOS = Orbs_nOrbsInState

  allocate (s(nOS, nOS))
  do j = 1, nOS
    do i = 1, nOS
      s(i, j) = 0.0_R64
      if (Method_Mb_OrbBased_bodyTypeOfOrb(i) .ne. Method_Mb_OrbBased_bodyTypeOfOrb(j)) cycle
      s(i, j) = Grid_InnerProduct(orbs(:, i), orbs(:, j))
    end do
  end do

  call LapackLib_DiagonalizeGeneric(evals, evecs, s, .true.)

  allocate (sInvSqrt(nOS, nOS))
  sInvSqrt = 0.0_R64
  do i = 1, nOS
    sInvSqrt(i, i) = 1.0_R64 / sqrt(evals(i))
  end do
  sInvSqrt = matmul(evecs, matmul(sInvSqrt, conjg(transpose(evecs))))

  call TransformOrbitalsKeepingState(coeffs, orbs, sInvSqrt)

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> FEDVR metric weights of a packed state [many-body part (nHead entries,
!> unit weights) | orbitals (Grid_Ylm_weights each)] for the SIL integrator.
function MetricWeights(nHead) result(res)
  use M_Grid
  use M_Grid_Ylm
  use M_Orbs

  integer(I32), intent(in) :: nHead
  real(R64), allocatable :: res(:)

  integer(I32) :: iOrb

  allocate (res(nHead + Grid_nPoints * Orbs_nOrbsInState))
  res(1:nHead) = 1.0_R64
  do iOrb = 1, Orbs_nOrbsInState
    res(nHead + (iOrb - 1) * Grid_nPoints + 1:nHead + iOrb * Grid_nPoints) = Grid_Ylm_weights(:)
  end do

end function
