! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Coeffs_Generic) S_Coeffs_Generic

  implicit none

!=============================================================================
! local procedures kinds
!=============================================================================

  complex(R64), allocatable, save :: coeffsBuffer(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Coeffs_Generic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Coeffs
    use M_ConfigList

    integer(I32) :: i

    call Say_Fabricate("coeffs.generic")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Coeffs_nCoeffs = 1
    do i = 1, size(configList)
      Coeffs_nCoeffs = Coeffs_nCoeffs * configList(i) % e % nConfigurations
    end do

    Coeffs_ApplyH1FillRdm1 => ApplyH1FillRdm1
    Coeffs_ApplyH2FillRdm2 => ApplyH2FillRdm2
    Coeffs_ApplyExcitation => ApplyExcitation
    Coeffs_IndexFromConfigurations => IndexFromConfigurations
    Coeffs_ConfigurationsFromIndex => ConfigurationsFromIndex

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyH1FillRdm1(coeffs, rdm1_, dCoeffs_, h1_)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Coeffs
    use M_ConfigList

    complex(R64), intent(in), contiguous, target                 :: coeffs(:)
    complex(R64), intent(out), allocatable, optional             :: rdm1_(:, :)
    complex(R64), intent(inout), contiguous, target, optional    :: dCoeffs_(:)
    complex(R64), intent(in), contiguous, optional               :: h1_(:, :)

    integer(I32) :: iBt1C, iCoeff, jC
    integer(I32) :: k1, iBt1, iTmp, orbCode, nO
    complex(R64) :: factor
    complex(R64), allocatable :: h1Lin(:)
    complex(R64), allocatable :: rdm1Lin(:)
    complex(R64), allocatable :: rdm1LinThread(:)
    integer(I32) :: configurations(Method_Mb_nBodyTypes)

    nO = Method_Mb_OrbBased_nOrbsSum

    if (present(rdm1_)) allocate (rdm1Lin(nO**2))
    if (present(rdm1_)) rdm1Lin(:) = 0.0_R64
    if (present(h1_)) allocate (h1Lin(nO**2))
    if (present(h1_)) call FillH1Lin(h1Lin, h1_)

    !$omp parallel default(shared) private(iCoeff, iTmp, iBt1, iBt1C, configurations, k1, orbCode, factor, jC, rdm1LinThread)

    if (present(rdm1_)) allocate (rdm1LinThread(nO**2))
    if (present(rdm1_)) rdm1LinThread(:) = 0.0_R64

    !$omp do
    do iCoeff = 1, Coeffs_nCoeffs

      iTmp = 1

      call Coeffs_ConfigurationsFromIndex(configurations, iCoeff)

      do iBt1 = 1, Method_Mb_nBodyTypes

        ! This multibase decomposition extracts the iBt1 component from the IndexFromConfigurationsd coeff array element iCoeff
        iBt1C = configurations(iBt1)

        do k1 = 1, configList(iBt1) % e % singles % nConnected(iBt1C)

          orbCode = configList(iBt1) % e % singles % orbCode(k1, iBt1C)
          factor = configList(iBt1) % e % singles % factor(k1, iBt1C)
          jC = configList(iBt1) % e % singles % excitedC(k1, iBt1C)

          ! This multibase transform turns the updated iBt1 component into the IndexFromConfigurationsd coeff array element jC
          jC = (jC - iBt1C)
          jC = iCoeff + jC * iTmp

          if (present(dCoeffs_) .and. present(h1_)) dCoeffs_(iCoeff) = dCoeffs_(iCoeff) + factor * h1Lin(orbCode) * coeffs(jC)
          if (present(rdm1_)) rdm1LinThread(orbCode) = rdm1LinThread(orbCode) + factor * conjg(coeffs(jC)) * coeffs(iCoeff)

        end do
        iTmp = iTmp * configList(iBt1) % e % nConfigurations

      end do
    end do
    !$omp end do

    !$omp critical
    if (present(rdm1_)) rdm1Lin(:) = rdm1Lin(:) + rdm1LinThread(:)
    !$omp end critical

    !$omp end parallel

    if (present(rdm1_)) then
      if (.not. allocated(rdm1_)) allocate (rdm1_(nO, nO))
      call FillRdm1FromRdm1Lin(rdm1_, rdm1Lin)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyH2FillRdm2(coeffs, rdm2_, dCoeffs_, h2_)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Coeffs
    use M_ConfigList

    complex(R64), intent(in), contiguous, target                 :: coeffs(:)
    complex(R64), intent(out), allocatable, optional             :: rdm2_(:, :, :, :)
    complex(R64), intent(inout), contiguous, target, optional    :: dCoeffs_(:)
    complex(R64), intent(in), contiguous, optional               :: h2_(:, :, :, :)

    integer(I32) :: iBt1C, iBt2C
    integer(I32) :: iCoeff, jC, jCTmp
    integer(I32) :: k1, k2, iBt1, iBt2, iTmp, jTmp, orbCode, orbCodeTmp
    integer(I32) :: nO
    complex(R64) :: factor, factorTmp

    complex(R64), allocatable :: h2Lin(:)
    complex(R64), allocatable :: rdm2Lin(:)
    complex(R64), allocatable :: rdm2LinThreads(:)
    integer(I32) :: configurations(Method_Mb_nBodyTypes)

    nO = Method_Mb_OrbBased_nOrbsSum

    if (present(rdm2_)) allocate (rdm2Lin(nO**4))
    if (present(rdm2_)) rdm2Lin(:) = 0.0_R64
    if (present(h2_)) allocate (h2Lin(nO**4))
    if (present(h2_)) call FillH2Lin(h2Lin, h2_)

    !$omp parallel default(shared) private(iCoeff, iTmp, iBt1, iBt1C, configurations, k1, orbCode, factor, jC, orbCodeTmp, &
    !$omp                                  factorTmp, jCTmp, jTmp, iBt2, iBt2C, k2, rdm2LinThreads)

    if (present(rdm2_)) allocate (rdm2LinThreads(nO**4))
    if (present(rdm2_)) rdm2LinThreads(:) = 0.0_R64

    !$omp do
    do iCoeff = 1, Coeffs_nCoeffs

      iTmp = 1

      call Coeffs_ConfigurationsFromIndex(configurations, iCoeff)

      do iBt1 = 1, Method_Mb_nBodyTypes

        ! This multibase decomposition extracts the iBt1 component from the IndexFromConfigurationsd coeff array element iCoeff
        iBt1C = configurations(iBt1)

        ! intra body type interaction
        do k1 = 1, configList(iBt1) % e % doubles % nConnected(iBt1C)

          orbCode = configList(iBt1) % e % doubles % orbCode(k1, iBt1C)
          factor = configList(iBt1) % e % doubles % factor(k1, iBt1C)
          jC = configList(iBt1) % e % doubles % excitedC(k1, iBt1C)

          ! This multibase transform turns the updated iBt1 component into the IndexFromConfigurationsd coeff array element jC
          jC = (jC - iBt1C)
          jC = iCoeff + jC * iTmp

          if (present(dCoeffs_)) dCoeffs_(iCoeff) = dCoeffs_(iCoeff) + factor * h2Lin(orbCode) * coeffs(jC)
          if (present(rdm2_)) rdm2LinThreads(orbCode) = rdm2LinThreads(orbCode) + factor * conjg(coeffs(jC)) * coeffs(iCoeff)

        end do

        if (iBt1 .eq. Method_Mb_nBodyTypes) cycle

        ! inter body type interaction
        do k1 = 1, configList(iBt1) % e % singles % nConnected(iBt1C)

          orbCodeTmp = configList(iBt1) % e % singles % orbCode(k1, iBt1C)
          orbCodeTmp = (orbCodeTmp - 1) * nO * nO

          factorTmp = configList(iBt1) % e % singles % factor(k1, iBt1C)
          jCTmp = configList(iBt1) % e % singles % excitedC(k1, iBt1C)

          ! This multibase transform turns the updated iBt1 component into the IndexFromConfigurationsd coeff array element jC
          jCTmp = (jCTmp - iBt1C)
          jCTmp = iCoeff + jCTmp * iTmp

          jTmp = iTmp * configList(iBt1) % e % nConfigurations

          do iBt2 = iBt1 + 1, Method_Mb_nBodyTypes

            ! This multibase decomposition extracts the iBt2 component from the IndexFromConfigurationsd coeff array element iCoeff
            iBt2C = configurations(iBt2)

            do k2 = 1, configList(iBt2) % e % singles % nConnected(iBt2C)

              orbCode = orbCodeTmp + configList(iBt2) % e % singles % orbCode(k2, iBt2C)
              factor = factorTmp * configList(iBt2) % e % singles % factor(k2, iBt2C)
              jC = configList(iBt2) % e % singles % excitedC(k2, iBt2C)

              ! This multibase transform turns the updated iBt1 component into the IndexFromConfigurationsd coeff array element jC
              jC = (jC - iBt2C)
              jC = jCTmp + jC * jTmp

              if (present(dCoeffs_)) dCoeffs_(iCoeff) = dCoeffs_(iCoeff) + factor * h2Lin(orbCode) * coeffs(jC)
              if (present(rdm2_)) rdm2LinThreads(orbCode) = rdm2LinThreads(orbCode) + factor * conjg(coeffs(jC)) * coeffs(iCoeff)

            end do

            jTmp = jTmp * configList(iBt2) % e % nConfigurations
          end do

        end do

        iTmp = iTmp * configList(iBt1) % e % nConfigurations
      end do

    end do
    !$omp end do

    !$omp critical
    if (present(rdm2_)) rdm2Lin(:) = rdm2Lin(:) + rdm2LinThreads(:)
    !$omp end critical

    !$omp end parallel

    if (present(rdm2_)) then
      if (.not. allocated(rdm2_)) allocate (rdm2_(nO, nO, nO, nO))
      call FillRdm2FromRdm2Lin(rdm2_, rdm2Lin)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm1FromRdm1Lin(rdm1_, rdm1Lin)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: rdm1_(:, :)
    complex(R64), intent(in), contiguous   :: rdm1Lin(:)

    integer(I32) :: i1, j1
    integer(I32) :: nO, orbCode

    nO = Method_Mb_OrbBased_nOrbsSum

    do j1 = 1, nO
      do i1 = 1, nO

        orbCode = 1
        orbCode = (orbCode - 1) * nO + i1
        orbCode = (orbCode - 1) * nO + j1

        rdm1_(i1, j1) = rdm1Lin(orbCode)

      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm2FromRdm2Lin(rdm2_, rdm2Lin)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: rdm2_(:, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm2Lin(:)

    integer(I32) :: i1, i2, j1, j2, i1bt, i2bt, j1bt, j2bt
    integer(I32) :: nO, orbCode

    nO = Method_Mb_OrbBased_nOrbsSum

    rdm2_(:, :, :, :) = 0.0_R64

    do j2 = 1, nO
      do j1 = 1, j2
        do i2 = 1, nO
          do i1 = 1, i2

            j2bt = Method_Mb_OrbBased_bodyTypeOfOrb(j2)
            j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)
            i2bt = Method_Mb_OrbBased_bodyTypeOfOrb(i2)
            i1bt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)

            if (i1bt .eq. i2bt .and. j1bt .ne. j2bt) cycle
            if (i1bt .ne. i2bt .and. j1bt .eq. j2bt) cycle

            orbCode = 1
            orbCode = (orbCode - 1) * nO + i1
            orbCode = (orbCode - 1) * nO + j1
            orbCode = (orbCode - 1) * nO + i2
            orbCode = (orbCode - 1) * nO + j2

            if (i1bt .eq. i2bt) then

              if (Method_Mb_bodyStatistics(i1bt) .eq. 'f') then

                rdm2_(i1, i2, j1, j2) = rdm2_(i1, i2, j1, j2) + rdm2Lin(orbCode)
                rdm2_(i2, i1, j1, j2) = rdm2_(i2, i1, j1, j2) - rdm2Lin(orbCode)
                rdm2_(i1, i2, j2, j1) = rdm2_(i1, i2, j2, j1) - rdm2Lin(orbCode)
                rdm2_(i2, i1, j2, j1) = rdm2_(i2, i1, j2, j1) + rdm2Lin(orbCode)

              else if (Method_Mb_bodyStatistics(i1bt) .eq. 'b') then

                rdm2_(i1, i2, j1, j2) = rdm2_(i1, i2, j1, j2) + rdm2Lin(orbCode)
                rdm2_(i2, i1, j1, j2) = rdm2_(i2, i1, j1, j2) + rdm2Lin(orbCode)
                rdm2_(i1, i2, j2, j1) = rdm2_(i1, i2, j2, j1) + rdm2Lin(orbCode)
                rdm2_(i2, i1, j2, j1) = rdm2_(i2, i1, j2, j1) + rdm2Lin(orbCode)

                if (i1 .eq. i2) rdm2_(i1, i2, j1, j2) = rdm2_(i1, i2, j1, j2) / 2.0_R64
                if (j1 .eq. j2) rdm2_(i1, i2, j1, j2) = rdm2_(i1, i2, j1, j2) / 2.0_R64

              else
                error stop "Method_Mb_bodyStatistics not implemented"
              end if

            else

              ! due to the klein transform we can choose the statistics of distinguishable
              ! here we choose antisymmetry
              rdm2_(i1, i2, j1, j2) = rdm2_(i1, i2, j1, j2) + rdm2Lin(orbCode)
              rdm2_(i2, i1, j1, j2) = rdm2_(i2, i1, j1, j2) - rdm2Lin(orbCode)
              rdm2_(i1, i2, j2, j1) = rdm2_(i1, i2, j2, j1) - rdm2Lin(orbCode)
              rdm2_(i2, i1, j2, j1) = rdm2_(i2, i1, j2, j1) + rdm2Lin(orbCode)

            end if

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH1Lin(h1Lin, h1_)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: h1Lin(:)
    complex(R64), intent(in), contiguous   :: h1_(:, :)

    integer(I32) :: i1, j1
    integer(I32) :: nO, orbCode

    nO = Method_Mb_OrbBased_nOrbsSum

    do j1 = 1, nO
      do i1 = 1, nO

        orbCode = 1
        orbCode = (orbCode - 1) * nO + i1
        orbCode = (orbCode - 1) * nO + j1

        h1Lin(orbCode) = h1_(i1, j1)

      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH2Lin(h2Lin, h2_)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: h2Lin(:)
    complex(R64), intent(in), contiguous   :: h2_(:, :, :, :)

    integer(I32) :: i1, i2, j1, j2, i1bt, i2bt, j1bt, j2bt
    integer(I32) :: nO, orbCode

    nO = Method_Mb_OrbBased_nOrbsSum

    do j2 = 1, nO
      do j1 = 1, j2
        do i2 = 1, nO
          do i1 = 1, i2

            j2bt = Method_Mb_OrbBased_bodyTypeOfOrb(j2)
            j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)
            i2bt = Method_Mb_OrbBased_bodyTypeOfOrb(i2)
            i1bt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)

            if (i1bt .eq. i2bt .and. j1bt .ne. j2bt) cycle
            if (i1bt .ne. i2bt .and. j1bt .eq. j2bt) cycle

            orbCode = 1
            orbCode = (orbCode - 1) * nO + i1
            orbCode = (orbCode - 1) * nO + j1
            orbCode = (orbCode - 1) * nO + i2
            orbCode = (orbCode - 1) * nO + j2

            if (i1bt .eq. i2bt) then

              if (Method_Mb_bodyStatistics(i1bt) .eq. 'f') then

                h2Lin(orbCode) = (h2_(i1, i2, j1, j2) - h2_(i1, i2, j2, j1) - h2_(i2, i1, j1, j2) + h2_(i2, i1, j2, j1))
                if (i1 .eq. i2) h2Lin(orbCode) = 0
                if (j1 .eq. j2) h2Lin(orbCode) = 0

              else if (Method_Mb_bodyStatistics(i1bt) .eq. 'b') then

                h2Lin(orbCode) = (h2_(i1, i2, j1, j2) + h2_(i1, i2, j2, j1) + h2_(i2, i1, j1, j2) + h2_(i2, i1, j2, j1))

                if (i1 .eq. i2) h2Lin(orbCode) = h2Lin(orbCode) / 2
                if (j1 .eq. j2) h2Lin(orbCode) = h2Lin(orbCode) / 2

              else
                error stop "Method_Mb_bodyStatistics not implemented"
              end if

            else

              ! due to the klein transform we can choose the statistics of distinguishable
              ! here we choose antisymmetry
              h2Lin(orbCode) = (h2_(i1, i2, j1, j2) - h2_(i1, i2, j2, j1) - h2_(i2, i1, j1, j2) + h2_(i2, i1, j2, j1))
              if (i1 .eq. i2) h2Lin(orbCode) = 0
              if (j1 .eq. j2) h2Lin(orbCode) = 0

            end if

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyExcitation(coeffs, creates, destroys, bt)
    use M_Coeffs
    use M_Method_Mb
    use M_ConfigList

    complex(R64), intent(inout), contiguous  :: coeffs(:)
    integer(I32), intent(in), contiguous  :: creates(:)
    integer(I32), intent(in), contiguous  :: destroys(:)
    integer(I32), intent(in) :: bt

    ! New variables for optimization
    integer(I32) :: strideBt, nConfigsBt, nUpper, stepSize
    integer(I32) :: c, cNew, k, ibt_loop
    integer(I32) :: offsetOld, offsetNew, loopIndex
    real(R64)  :: factor

    ! Reuse buffer to avoid reallocation
    if (.not. allocated(coeffsBuffer)) allocate (coeffsBuffer(size(coeffs)))
    if (size(coeffsBuffer) .ne. size(coeffs)) then
      deallocate (coeffsBuffer)
      allocate (coeffsBuffer(size(coeffs)))
    end if

    coeffsBuffer = coeffs

    coeffs(:) = 0.0_R64

    ! Compute strides for tensor product structure
    strideBt = 1
    do ibt_loop = 1, bt - 1
      strideBt = strideBt * configList(ibt_loop) % e % nConfigurations
    end do

    nConfigsBt = configList(bt) % e % nConfigurations
    stepSize = strideBt * nConfigsBt
    nUpper = Coeffs_nCoeffs / stepSize

    ! Iterate only over the relevant body type configurations
    do c = 1, nConfigsBt

      ! Excite only the configuration for body type bt
      call configList(bt) % e % ExciteConfiguration(cNew, factor, creates, destroys, c)
      if (cNew .eq. 0) cycle

      offsetOld = (c - 1) * strideBt
      offsetNew = (cNew - 1) * strideBt

      ! Apply to all tensor blocks (vectorized copy)
      do k = 0, nUpper - 1
        loopIndex = k * stepSize
        coeffs(loopIndex + offsetNew + 1:loopIndex + offsetNew + strideBt) = &
          factor * coeffsBuffer(loopIndex + offsetOld + 1:loopIndex + offsetOld + strideBt)
      end do

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pure subroutine ConfigurationsFromIndex(configurations, iCoeff)
    use M_ConfigList

    integer(I32), intent(out), contiguous :: configurations(:)
    integer(I32), intent(in) :: iCoeff

    integer(I32) :: ibt, iTmp

    iTmp = iCoeff - 1
    do ibt = 1, size(configurations)

      ! This multibase decomposition extracts the ibt component from the IndexFromConfigurationsd coeff array element iCoeff
      configurations(ibt) = mod(iTmp, configList(ibt) % e % nConfigurations) + 1
      iTmp = iTmp / configList(ibt) % e % nConfigurations

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pure subroutine IndexFromConfigurations(iCoeff, configurations)
    use M_ConfigList

    integer(I32), intent(out) :: iCoeff
    integer(I32), intent(in), contiguous  :: configurations(:)

    integer(I32) :: ibt, iTmp

    iTmp = 1
    iCoeff = 0
    do ibt = 1, size(configurations)
      iCoeff = iCoeff + (configurations(ibt) - 1) * iTmp
      iTmp = iTmp * configList(ibt) % e % nConfigurations
    end do

    iCoeff = iCoeff + 1

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
