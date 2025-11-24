! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Coeffs_Hubbard_NoSpinSym) S_Coeffs_Hubbard_NoSpinSym

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Coeffs_Hubbard_NoSpinSym_Fabricate
    use M_Utils_Say
    use M_Coeffs
    use M_Coeffs_Hubbard

    call Say_Fabricate("coeffs.hubbard.noSpinSym")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Coeffs_ApplyH1FillRdm1 => ApplyH1FillRdm1
    Coeffs_ApplyH2FillRdm2 => ApplyH2FillRdm2
    Coeffs_ConfigurationsFromIndex => ConfigurationsFromIndex
    Coeffs_IndexFromConfigurations => IndexFromConfigurations

    Coeffs_nCoeffs = Coeffs_Hubbard_nCoeffsUP * Coeffs_Hubbard_nCoeffsDN

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyH1FillRdm1(coeffs, rdm1_, dCoeffs_, h1_)
    use M_Utils_UnusedVariables
    use M_Coeffs_Hubbard

    complex(R64), intent(in), contiguous, target                 :: coeffs(:)
    complex(R64), intent(out), allocatable, optional             :: rdm1_(:, :)
    complex(R64), intent(inout), contiguous, target, optional    :: dCoeffs_(:)
    complex(R64), intent(in), contiguous, optional               :: h1_(:, :)

    complex(R64), contiguous, pointer :: coeffsReshape(:, :) ! Probably "contiguous" is not necessary
    complex(R64), contiguous, pointer :: dCoeffsReshape(:, :)

    integer(I32) :: i1, i2, k, it
    real(R64) :: weight

    if (.false.) call UnusedVariables_Mark(h1_)
    if (present(rdm1_)) error stop "calculating rdm1_ in ApplyH1FillRdm1 not supported"

    coeffsReshape(1:Coeffs_Hubbard_nCoeffsUP, 1:Coeffs_Hubbard_nCoeffsDN) => coeffs(:)
    dCoeffsReshape(1:Coeffs_Hubbard_nCoeffsUP, 1:Coeffs_Hubbard_nCoeffsDN) => dCoeffs_(:)

    !==========================================================================================
    ! Kinetic part of the hamiltonian
    !==========================================================================================

    !$omp parallel default(shared) private(i1, i2, k, weight, it)

    !$omp do collapse(2)
    Do i2 = 1, Coeffs_Hubbard_nCoeffsDN
      Do i1 = 1, Coeffs_Hubbard_nCoeffsUP

        dCoeffsReshape(i1, i2) = 0.0_R64

        Do k = 1, Coeffs_Hubbard_nConnectedUP(i1)

          weight = Coeffs_Hubbard_weightUP(k, i1)
          it = Coeffs_Hubbard_hoppUP(k, i1)

          dCoeffsReshape(i1, i2) = dCoeffsReshape(i1, i2) + weight * coeffsReshape(it, i2)
        end do

        Do k = 1, Coeffs_Hubbard_NConnectedDN(i2)

          weight = Coeffs_Hubbard_weightDN(k, i2)
          it = Coeffs_Hubbard_hoppDN(k, i2)

          dCoeffsReshape(i1, i2) = dCoeffsReshape(i1, i2) + weight * coeffsReshape(i1, it)
        end do

      end do
    end do
    !$omp end do

    !$omp end parallel

    !==========================================================================================
    ! Potential part of the hamiltonian
    !==========================================================================================

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

    complex(R64), contiguous, pointer :: coeffsReshape(:, :) ! Probably "contiguous" is not necessary
    complex(R64), contiguous, pointer :: dCoeffsReshape(:, :)

    integer(I32) :: i1, i2, i
    integer(I64) :: iTmp
    real(R64) :: weight
    real(R64) :: hubU

    if (present(rdm2_)) error stop "calculating rdm2_ in hubbardCoeffs not supported"

    hubU = real(h2_(1, 1 + Method_Mb_OrbBased_nOrbs(1), 1, 1 + Method_Mb_OrbBased_nOrbs(1)), kind=R64) * 2.0_R64

    coeffsReshape(1:Coeffs_Hubbard_nCoeffsUP, 1:Coeffs_Hubbard_nCoeffsDN) => coeffs(:)
    dCoeffsReshape(1:Coeffs_Hubbard_nCoeffsUP, 1:Coeffs_Hubbard_nCoeffsDN) => dCoeffs_(:)

    !==========================================================================================
    ! Interaction part of the hamiltonian
    !==========================================================================================

    !$omp parallel default(shared) private(i1, i2, iTmp, weight, i)

    !$omp do collapse(2)
    Do i2 = 1, Coeffs_Hubbard_nCoeffsDN
      Do i1 = 1, Coeffs_Hubbard_nCoeffsUP

        !-------------------------------------------------------------------------------------
        ! calculate number of doubly occupied Grid_nPoints for each configuration:
        ! ovlap = popcnt(iand(Coeffs_Hubbard_bitcodesUP(i1), Coeffs_Hubbard_bitcodesDN(i2)))
        ! is the number of doubly occupied Grid_nPoints for combining configuration i1 and i2
        !-------------------------------------------------------------------------------------

        iTmp = iand(Coeffs_Hubbard_bitcodesUP(i1), Coeffs_Hubbard_bitcodesDN(i2))

        weight = 0.0_R64
        do i = 1, Grid_nPoints
          if (btest(iTmp, i - 1)) weight = weight + hubU
        end do

        dCoeffsReshape(i1, i2) = dCoeffsReshape(i1, i2) + weight * coeffsReshape(i1, i2)

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

    configurations(1) = Mod(iCoeff - 1, Coeffs_Hubbard_nCoeffsUP) + 1
    configurations(2) = (iCoeff - 1) / Coeffs_Hubbard_nCoeffsUP + 1

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pure subroutine IndexFromConfigurations(iCoeff, configurations)
    use M_Coeffs_Hubbard
    use M_Method

    integer(I32), intent(out) :: iCoeff
    integer(I32), intent(in), contiguous  :: configurations(:)

    iCoeff = (configurations(2) - 1) * Coeffs_Hubbard_nCoeffsUP + configurations(1) ! -1 +1

  end subroutine

end submodule
