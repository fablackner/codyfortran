! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Basis change of ket-first RDMs, rdm2(i1,i2,j1,j2) = <a+_j1 a+_j2 a_i2 a_i1>
!> (rank 6 alike).
module M_Utils_RdmBasisChange
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> rdm <- m rdm m^+ in every particle: x(.., i, ..) <- sum_k m(i, k)
  !> x(.., k, ..) on the r/2 ket indices, conjg(m) on the r/2 bra indices of a
  !> ket-first RDM of even rank r and extent n = size(m, 1), column-major in x.
  !> m = U^+ gives the RDM in the orbitals phi U, m = <phi(t)|phi_ref> moves a
  !> kernel from the reference orbitals to the current ones, m = exp(-i k z)
  !> applies the one-body unitary exp(-i k Z) to the state. One ZGEMM per
  !> index (threaded BLAS: call outside parallel regions): it transforms the
  !> leading index and moves it last, so r of them restore the order.
  subroutine RdmBasisChange_Transform(x, r, m)
    use M_Utils_BlasLib
    use M_Utils_Constants, only: CONE, CZERO

    complex(R64), intent(inout)          :: x(*)
    integer(I32), intent(in)             :: r
    complex(R64), intent(in), contiguous :: m(:, :)

    complex(R64), allocatable :: y(:), mBra(:, :)
    integer(I32) :: n, mode

    if (mod(r, 2) .ne. 0) error stop "RdmBasisChange_Transform: the rank must be even"
    n = size(m, 1)
    mBra = conjg(m)
    allocate (y(n**r))
    do mode = 1, r, 2
      call Step(y, x, mode)
      call Step(x, y, mode + 1)
    end do

  contains

    !> dst(rest, i) = sum_k a(i, k) src(k, rest), a = m on the ket and
    !> conjg(m) on the bra indices
    subroutine Step(dst, src, mode)
      complex(R64), intent(inout) :: dst(*)
      complex(R64), intent(in)    :: src(*)
      integer(I32), intent(in)    :: mode

      if (mode <= r / 2) then
        call ZGEMM('T', 'T', n**(r - 1), n, n, CONE, src, n, m, n, CZERO, dst, n**(r - 1))
      else
        call ZGEMM('T', 'T', n**(r - 1), n, n, CONE, src, n, mBra, n, CZERO, dst, n**(r - 1))
      end if

    end subroutine

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> A block of a ket-first RDM of even rank r, given in a subset of orbitals,
  !> in all n orbitals: y(.., i, ..) = sum_k mKet(i, k) x(.., k, ..) on the r/2
  !> ket indices (extent size(mKet, 2) in x) and conjg(mBra) on the r/2 bra
  !> indices (extent size(mBra, 2)); y has n^r elements, n = size(mKet, 1).
  !> With the columns of two sets of orbitals as mKet and mBra this embeds
  !> the block that annihilates in the first and creates in the second set.
  !> One ZGEMM per index, as RdmBasisChange_Transform, on arrays that grow
  !> from the block to n^r.
  subroutine RdmBasisChange_EmbedBlock(y, x, r, mKet, mBra)
    use M_Utils_BlasLib
    use M_Utils_Constants, only: CONE, CZERO

    complex(R64), intent(out)            :: y(*)
    complex(R64), intent(in)             :: x(*)
    integer(I32), intent(in)             :: r
    complex(R64), intent(in), contiguous :: mKet(:, :)
    complex(R64), intent(in), contiguous :: mBra(:, :)

    complex(R64), allocatable :: work(:), mBraConjg(:, :)
    integer(I32) :: n, nKet, nBra, mode, nRest

    if (mod(r, 2) .ne. 0) error stop "RdmBasisChange_EmbedBlock: the rank must be even"
    n = size(mKet, 1)
    nKet = size(mKet, 2)
    nBra = size(mBra, 2)
    mBraConjg = conjg(mBra)
    allocate (work(n**r))

    ! dst(rest, i) = sum_k a(i, k) src(k, rest): the leading index is
    ! transformed and moved last; r steps (even) end in y
    do mode = 1, r
      if (mode <= r / 2) then
        nRest = nKet**(r / 2 - mode) * nBra**(r / 2) * n**(mode - 1)
        if (mode .eq. 1) then
          call ZGEMM('T', 'T', nRest, n, nKet, CONE, x, nKet, mKet, n, CZERO, work, nRest)
        else if (mod(mode, 2) .eq. 0) then
          call ZGEMM('T', 'T', nRest, n, nKet, CONE, work, nKet, mKet, n, CZERO, y, nRest)
        else
          call ZGEMM('T', 'T', nRest, n, nKet, CONE, y, nKet, mKet, n, CZERO, work, nRest)
        end if
      else
        nRest = nBra**(r - mode) * n**(mode - 1)
        if (mod(mode, 2) .eq. 0) then
          call ZGEMM('T', 'T', nRest, n, nBra, CONE, work, nBra, mBraConjg, n, CZERO, y, nRest)
        else
          call ZGEMM('T', 'T', nRest, n, nBra, CONE, y, nBra, mBraConjg, n, CZERO, work, nRest)
        end if
      end if
    end do

  end subroutine

end module
