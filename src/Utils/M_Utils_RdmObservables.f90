! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Core routines for assembling and manipulating reduced density matrices.
!>
!> Provides 1/2/3-body RDM construction, projections, contractions, and
!> normalization/trace utilities used throughout the analysis pipeline.
module M_Utils_RdmObservables
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function RdmObservables_Energy(rdm1, h1, rdm2, h2) result(res)

    real(R64)                :: res
    complex(R64), intent(in), contiguous  :: rdm1(:, :)
    complex(R64), intent(in), contiguous  :: h1(:, :)
    complex(R64), intent(in), contiguous  :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)

    complex(R64) :: cTmp

    integer(I32) :: k1, l1, k2, l2, nO

    nO = size(h1, 1)

    cTmp = 0.0_R64
    do k1 = 1, nO
      do l1 = 1, nO
        cTmp = cTmp + h1(k1, l1) * rdm1(l1, k1)
      end do
    end do

    do k1 = 1, nO
      do k2 = 1, nO
        do l1 = 1, nO
          do l2 = 1, nO

            cTmp = cTmp + h2(k1, k2, l1, l2) * rdm2(l1, l2, k1, k2)

          end do
        end do
      end do
    end do

    res = real(cTmp, kind=R64)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function RdmObservables_Norm(rdm2) result(norm)
    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    real(R64) :: norm
    integer(I32) :: nO, i1, i2
    complex(R64) :: cNorm

    nO = size(rdm2, 1)
    cNorm = 0.0_R64
    do i1 = 1, nO
      do i2 = 1, nO
        cNorm = cNorm + rdm2(i1, i2, i1, i2)
      end do
    end do
    norm = real(cNorm, kind=R64)
  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function RdmObservables_Spin(rdm2) result(spin)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    real(R64) :: spin
    integer(I32) :: nOa, i1, i2, nBodiesUP
    complex(R64) :: cSpin

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nBodiesUP = Method_Mb_nBodies(1)

    cSpin = 0.0_R64
    do i1 = 1, nOa
      do i2 = 1, nOa
        cSpin = cSpin + rdm2(i1, i2 + nOa, i2, i1 + nOa)
      end do
    end do
    spin = real(nBodiesUP - cSpin, kind=R64)
  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function RdmObservables_SingletTraceCondition(rdm2) result(spinTr1)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    real(R64) :: spinTr1
    integer(I32) :: nOa, i1, j1, i2
    complex(R64) :: norm1, norm2

    nOa = Method_Mb_OrbBased_nOrbs(1)

    spinTr1 = 0.0_R64
    do i1 = 1, nOa
      do j1 = 1, nOa
        norm1 = 0.0_R64
        norm2 = 0.0_R64
        do i2 = 1, nOa
          norm1 = norm1 + rdm2(i1, i2 + nOa, i1, i2 + nOa)
          norm2 = norm2 + rdm2(i1, i2 + nOa, i2, i1 + nOa)
        end do
        spinTr1 = spinTr1 + abs(norm1 - Method_Mb_nBodies(1) * norm2)
      end do
    end do
  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function RdmObservables_PermutationSymmetry(rdm2) result(permSym)
    use M_Method_Mb_OrbBased
    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    real(R64) :: permSym
    integer(I32) :: nO, nOa, i1, i2, j1, j2

    nO = size(rdm2, 1)
    nOa = Method_Mb_OrbBased_nOrbs(1)

    permSym = 0.0_R64
    do j2 = 1, nO
      do j1 = 1, nO
        do i2 = 1, nO
          do i1 = 1, nO
            permSym = permSym + &
                      abs(rdm2(i1, i2, j1, j2) + rdm2(i2, i1, j1, j2)) + &
                      abs(rdm2(i1, i2, j1, j2) + rdm2(i1, i2, j2, j1)) + &
                      abs(rdm2(i1, i2, j1, j2) - rdm2(i2, i1, j2, j1))
          end do
        end do
      end do
    end do
    permSym = permSym / (real(nOa, kind=R64))**4
  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function RdmObservables_Hermiticity(rdm2) result(hermSym)
    use M_Method_Mb_OrbBased
    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    real(R64) :: hermSym
    integer(I32) :: nO, nOa, i1, i2, j1, j2

    nO = size(rdm2, 1)
    nOa = Method_Mb_OrbBased_nOrbs(1)

    hermSym = 0.0_R64
    do j2 = 1, nO
      do j1 = 1, nO
        do i2 = 1, nO
          do i1 = 1, nO
            hermSym = hermSym + abs(rdm2(i1, i2, j1, j2) - conjg(rdm2(j1, j2, i1, i2)))
          end do
        end do
      end do
    end do
    hermSym = hermSym / (real(nOa, kind=R64))**4
  end function

end module
