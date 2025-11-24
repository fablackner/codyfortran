! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Coeffs_Hubbard_PlusSpinSym) S_Coeffs_Hubbard_PlusSpinSym

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Coeffs_Hubbard_PlusSpinSym_Fabricate
    use M_Utils_Say
    use M_Coeffs
    use M_Coeffs_Hubbard

    call Say_Fabricate("coeffs.hubbard.plusSpinSym")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Coeffs_ApplyH1FillRdm1 => ApplyH1FillRdm1
    Coeffs_ApplyH2FillRdm2 => ApplyH2FillRdm2
    Coeffs_ConfigurationsFromIndex => ConfigurationsFromIndex
    Coeffs_IndexFromConfigurations => IndexFromConfigurations

    Coeffs_nCoeffs = Coeffs_Hubbard_nCoeffsUP * (Coeffs_Hubbard_nCoeffsUP + 1) / 2

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyH1FillRdm1(coeffs, rdm1_, dCoeffs_, h1_)
    use M_Utils_UnusedVariables
    use M_Coeffs_Hubbard

    complex(R64), intent(in), contiguous, target                 :: coeffs(:)
    complex(R64), intent(out), allocatable, optional             :: rdm1_(:, :)
    complex(R64), intent(inout), contiguous, target, optional    :: dCoeffs_(:)
    complex(R64), intent(in), contiguous, optional               :: h1_(:, :)

    complex(R64) :: coeffsReshape(1:Coeffs_Hubbard_nCoeffsUP, 1:Coeffs_Hubbard_nCoeffsUP) ! Probably "contiguous" is not necessary
    complex(R64) :: dCoeffsReshape(1:Coeffs_Hubbard_nCoeffsUP, 1:Coeffs_Hubbard_nCoeffsUP)

    integer(I32) :: i1, i2, k, it, iCoeff
    real(R64) :: weight

    if (.false.) call UnusedVariables_Mark(h1_)
    if (present(rdm1_)) error stop "calculating rdm1_ in ApplyH1FillRdm1 not supported"

    !$omp parallel default(shared) private(i1, i2, iCoeff)

    !$omp do
    Do i2 = 1, Coeffs_Hubbard_nCoeffsUP
      Do i1 = 1, i2 - 1
        iCoeff = (i2 - 1) * i2 / 2 + i1

        coeffsReshape(i1, i2) = coeffs(iCoeff) / sqrt(2.0_R64)

      end do
      iCoeff = (i2 - 1) * i2 / 2 + i2

      coeffsReshape(i2, i2) = coeffs(iCoeff)

    end do
    !$omp end do

    !$omp end parallel

    !$omp parallel default(shared) private(i1, i2)

    !$omp do
    Do i2 = 1, Coeffs_Hubbard_nCoeffsUP
      Do i1 = 1, i2 - 1
        coeffsReshape(i2, i1) = coeffsReshape(i1, i2)
      end do
    end do
    !$omp end do

    !$omp end parallel

    !==========================================================================================
    ! Kinetic part of the hamiltonian
    !==========================================================================================

    !$omp parallel default(shared) private(i2, k, weight, it)

    !$omp do
    Do i2 = 1, Coeffs_Hubbard_nCoeffsUP

      dCoeffsReshape(:, i2) = 0.0_R64

      Do k = 1, Coeffs_Hubbard_nConnectedUP(i2)

        weight = Coeffs_Hubbard_weightDN(k, i2)
        it = Coeffs_Hubbard_hoppDN(k, i2)

        dCoeffsReshape(:, i2) = dCoeffsReshape(:, i2) + weight * coeffsReshape(:, it)

      end do

    end do
    !$omp end do

    !$omp end parallel

    !$omp parallel default(shared) private(i1, i2)

    !$omp do
    Do i2 = 1, Coeffs_Hubbard_nCoeffsUP
      Do i1 = 1, i2
        dCoeffsReshape(i1, i2) = dCoeffsReshape(i1, i2) + dCoeffsReshape(i2, i1)
        dCoeffsReshape(i2, i1) = dCoeffsReshape(i1, i2)
      end do
    end do
    !$omp end do

    !$omp end parallel

    !$omp parallel default(shared) private(i1, i2, iCoeff)

    !$omp do
    Do i2 = 1, Coeffs_Hubbard_nCoeffsUP
      Do i1 = 1, i2 - 1
        iCoeff = (i2 - 1) * i2 / 2 + i1

        dCoeffs_(iCoeff) = dCoeffsReshape(i1, i2) * sqrt(2.0_R64)

      end do
      iCoeff = (i2 - 1) * i2 / 2 + i2

      dCoeffs_(iCoeff) = dCoeffsReshape(i1, i2)

    end do
    !$omp end do

    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyH2FillRdm2(coeffs, rdm2_, dCoeffs_, h2_)
    use M_Grid
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Coeffs_Hubbard

    complex(R64), intent(in), contiguous, target                 :: coeffs(:)
    complex(R64), intent(out), allocatable, optional             :: rdm2_(:, :, :, :)
    complex(R64), intent(inout), contiguous, target, optional    :: dCoeffs_(:)
    complex(R64), intent(in), contiguous, optional               :: h2_(:, :, :, :)

    integer(I32) :: i1, i2, iCoeff, i
    integer(I64) :: iTmp
    real(R64) :: weight
    real(R64) :: hubU

    if (present(rdm2_)) error stop "calculating rdm2_ in hubbardCoeffs not supported"

    hubU = real(h2_(1, 1 + Method_Mb_OrbBased_nOrbs(1), 1, 1 + Method_Mb_OrbBased_nOrbs(1)), kind=R64) * 2.0_R64

    !==========================================================================================
    ! Interaction part of the hamiltonian
    !==========================================================================================

    !$omp parallel default(shared) private(i1, i2, iCoeff, iTmp, weight, i)

    !$omp do
    Do i2 = 1, Coeffs_Hubbard_nCoeffsUP
      Do i1 = 1, i2
        iCoeff = (i2 - 1) * i2 / 2 + i1

        iTmp = iand(Coeffs_Hubbard_bitcodesUP(i1), Coeffs_Hubbard_bitcodesUP(i2))

        weight = 0.0_R64
        do i = 1, Grid_nPoints
          if (btest(iTmp, i - 1)) weight = weight + hubU
        end do

        dCoeffs_(iCoeff) = dCoeffs_(iCoeff) + weight * coeffs(iCoeff)

      end do

    end do
    !$omp end do

    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pure subroutine ConfigurationsFromIndex(configurations, iCoeff)
    use M_Coeffs_Hubbard
    use M_Method

    integer(I32), intent(out), contiguous :: configurations(:)
    integer(I32), intent(in) :: iCoeff

    integer(I32) :: i2, tri

    Do i2 = Coeffs_Hubbard_nCoeffsUP, 1, -1
      tri = (i2) * (i2 + 1) / 2
      if (tri > iCoeff) cycle
      configurations(2) = i2
      configurations(1) = iCoeff - tri
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pure subroutine IndexFromConfigurations(iCoeff, configurations)
    use M_Coeffs_Hubbard
    use M_Method

    integer(I32), intent(out) :: iCoeff
    integer(I32), intent(in), contiguous  :: configurations(:)
    integer(I32) :: minVal, maxVal

    minVal = min(configurations(1), configurations(2))
    maxVal = max(configurations(1), configurations(2))

    iCoeff = (maxVal) * (maxVal + 1) / 2 + minVal

  end subroutine

end submodule
