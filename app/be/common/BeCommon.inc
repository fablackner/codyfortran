!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Helpers shared by the beryllium programs (included after `contains`).
!
! - EmbedOrbitals:        load orbXX_YY.in written on a grid with a smaller
!                         lmax (same radial grid and mmax)
! - HermitianFunction:    U diag(f) U^H
! - MetricWeights:        FEDVR metric of a packed state (SIL integrator)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fill the orbitals from orbXX_YY.in files written on a Ylm grid with the same
!> radial grid and mmax but a smaller lmax (lmaxIn). The channel layout is
!> l-major, so the stored orbital is a prefix of the current one and the
!> missing high-l channels are zero. Like OrbsInit_Initialize, it ends with
!> the projection on the declared symmetry (Orbs_ProjectOnSymmetry), which
!> removes round-off admixtures of other m from the loaded orbitals.
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

  if (lmaxIn > Grid_Ylm_lmax) error stop "EmbedOrbitals: lmaxIn exceeds grid.ylm.lmax"

  nIn = Grid_Ylm_NChannels(lmaxIn, Grid_Ylm_mmax) * Grid_Ylm_nRadial
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

  call Orbs_ProjectOnSymmetry(orbs)

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The Hermitian matrix U diag(f) U^H of eigenvectors U and values f.
function HermitianFunction(evecs, f) result(res)
  complex(R64), intent(in), contiguous :: evecs(:, :)
  real(R64), intent(in), contiguous    :: f(:)
  complex(R64), allocatable :: res(:, :)

  res = matmul(evecs, spread(f, 2, size(f)) * conjg(transpose(evecs)))

end function

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
