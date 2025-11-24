! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Coeffs_Hubbard) S_Coeffs_Hubbard

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  integer(I32), allocatable :: ConfigurationsFromIndexUP(:)
  integer(I32), allocatable :: ConfigurationsFromIndexDN(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Coeffs_Hubbard_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_SfGslLib
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Coeffs
    use M_Coeffs_Hubbard_NoSpinSym
    use M_Coeffs_Hubbard_PlusSpinSym
    use M_Coeffs_Hubbard_MinusSpinSym

    call Say_Fabricate("coeffs.hubbard")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Coeffs_Hubbard_nCoeffsUP = SfGslLib_Binomial(Method_Mb_OrbBased_nOrbs(1), Method_Mb_nBodies(1))
    Coeffs_Hubbard_nCoeffsDN = SfGslLib_Binomial(Method_Mb_OrbBased_nOrbs(2), Method_Mb_nBodies(2))

    Coeffs_Setup => Setup
    Coeffs_ApplyExcitation => ApplyExcitation

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("coeffs.hubbard.noSpinSym")) then
      call Coeffs_Hubbard_NoSpinSym_Fabricate

    else if (Json_GetExistence("coeffs.hubbard.plusSpinSym")) then
      call Coeffs_Hubbard_PlusSpinSym_Fabricate

    else if (Json_GetExistence("coeffs.hubbard.minusSpinSym")) then
      call Coeffs_Hubbard_MinusSpinSym_Fabricate

    else
      error stop "coeffs.hubbard is missing one of: noSpinSym, plusSpinSym, minusSpinSym"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say

    call Say_Setup("coeffs.hubbard")

    call SetupConfigurations
    call SetupHamiltonian

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine SetupConfigurations
    use M_Utils_CombinationGslLib
    use M_Grid
    use M_Method_Mb

    integer(I32) :: iCoeff, j
    integer(I32), allocatable :: collect(:, :)

    allocate (Coeffs_Hubbard_bitcodesUP(Coeffs_Hubbard_nCoeffsUP))
    allocate (Coeffs_Hubbard_bitcodesDN(Coeffs_Hubbard_nCoeffsDN))

    allocate (ConfigurationsFromIndexUP(0:2**(Grid_nPoints)))
    allocate (ConfigurationsFromIndexDN(0:2**(Grid_nPoints)))

    allocate (collect(Method_Mb_nBodies(1), Coeffs_Hubbard_nCoeffsUP))

    !-------------------------------------------------------------------------------------
    ! collect(i, iCoeff) is the position of the i-th particle in configuration iCoeff
    !-------------------------------------------------------------------------------------

    call CombinationGslLib_CombiNoRepeat(collect, Grid_nPoints)

    !-------------------------------------------------------------------------------------
    ! Coeffs_Hubbard_bitcodesUP(iCoeff) IndexFromConfigurationss the distribution of particles in configuration iCoeff
    ! the binary digit of Coeffs_Hubbard_bitcodesUP(iCoeff) is equal 1 for particles and 0 for holes
    !-------------------------------------------------------------------------------------

    Coeffs_Hubbard_bitcodesUP(:) = 0
    do j = 1, Method_Mb_nBodies(1)
      do iCoeff = 1, Coeffs_Hubbard_nCoeffsUP
        Coeffs_Hubbard_bitcodesUP(iCoeff) = ibset(Coeffs_Hubbard_bitcodesUP(iCoeff), collect(j, iCoeff) - 1)
      end do
    end do

    !-------------------------------------------------------------------------------------
    ! ConfigurationsFromIndexUP allows to calculate the configuration index iCoeff associated
    ! with a particluar distribution    ...see function inverse
    !-------------------------------------------------------------------------------------

    ConfigurationsFromIndexUP(:) = 0
    do iCoeff = 1, Coeffs_Hubbard_nCoeffsUP
      ConfigurationsFromIndexUP(Coeffs_Hubbard_bitcodesUP(iCoeff)) = iCoeff
    end do

    deallocate (collect)

    !-------------------------------------------------------------------------------------
    ! repeat for spin down
    !-------------------------------------------------------------------------------------

    allocate (collect(Method_Mb_nBodies(2), Coeffs_Hubbard_nCoeffsDN))
    if (Method_Mb_nBodies(2) .ne. 0) call CombinationGslLib_CombiNoRepeat(collect, Grid_nPoints)

    Coeffs_Hubbard_bitcodesDN(:) = 0
    do j = 1, Method_Mb_nBodies(2)
      do iCoeff = 1, Coeffs_Hubbard_nCoeffsDN
        Coeffs_Hubbard_bitcodesDN(iCoeff) = ibset(Coeffs_Hubbard_bitcodesDN(iCoeff), collect(j, iCoeff) - 1)
      end do
    end do

    ConfigurationsFromIndexDN(:) = 0
    do iCoeff = 1, Coeffs_Hubbard_nCoeffsDN
      ConfigurationsFromIndexDN(Coeffs_Hubbard_bitcodesDN(iCoeff)) = iCoeff
    end do

    deallocate (collect)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine SetupHamiltonian
    use M_Grid
    use M_Grid_Lattice
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Utils_Combinatorics

    integer(I32) :: iCoeff, i1, j1, count, k, iDn, iUp, iSite, bitpos
    integer(I64) :: bits, tbits, interactionBits
    integer(I32) :: nDim

    integer(I32), allocatable :: permut(:)
    real(R64), allocatable :: upTmp(:)
    real(R64), allocatable :: dnTmp(:)
    real(R64) :: interactionStrength, onSite

    ! Make sure Method_Mb_OrbBased_h1 is filled
    if (.not. allocated(Method_Mb_OrbBased_h1)) then
      error stop "Method_Mb_OrbBased_h1 is not allocated, please call Method_Mb_OrbBased_FillH1 first"
    end if

    ! Make sure Method_Mb_OrbBased_h2 is filled
    if (.not. allocated(Method_Mb_OrbBased_h2)) then
      error stop "Method_Mb_OrbBased_h2 is not allocated, please call Method_Mb_OrbBased_FillH2 first"
    end if

    allocate (Coeffs_Hubbard_nConnectedUP(Coeffs_Hubbard_nCoeffsUP))
    allocate (Coeffs_Hubbard_NConnectedDN(Coeffs_Hubbard_nCoeffsDN))

    nDim = 0
    if (Grid_Lattice_xSize > 1) nDim = nDim + 1
    if (Grid_Lattice_ySize > 1) nDim = nDim + 1
    if (Grid_Lattice_zSize > 1) nDim = nDim + 1

    allocate (Coeffs_Hubbard_hoppUP(Method_Mb_nBodies(1) * nDim * 2, Coeffs_Hubbard_nCoeffsUP))
    allocate (Coeffs_Hubbard_hoppDN(Method_Mb_nBodies(2) * nDim * 2, Coeffs_Hubbard_nCoeffsDN))

    allocate (Coeffs_Hubbard_weightUP(Method_Mb_nBodies(1) * nDim * 2, Coeffs_Hubbard_nCoeffsUP))
    allocate (Coeffs_Hubbard_weightDN(Method_Mb_nBodies(2) * nDim * 2, Coeffs_Hubbard_nCoeffsDN))

    ! Allocate arrays for onsite interaction
    allocate (Coeffs_Hubbard_interactionValues(2**Grid_nPoints))

    allocate (upTmp(Method_Mb_nBodies(1) * nDim * 2))
    allocate (dnTmp(Method_Mb_nBodies(1) * nDim * 2))

    !-------------------------------------------------------------------------------------
    ! calculate the configurations that are connected by hopping of one of the particles
    ! this allows for an efficient calculation of the kinetic part in the hamiltonian
    !-------------------------------------------------------------------------------------

    Coeffs_Hubbard_nConnectedUP(:) = 0
    Coeffs_Hubbard_nConnectedDN(:) = 0

    ! This assumes the orbs are in position states!
    Do i1 = 1, Grid_nPoints
      Do j1 = 1, Grid_nPoints
        ! Skip if there's no hopping between these sites.
        if (abs(Method_Mb_OrbBased_h1(i1, j1)) < 1e-15_R64) cycle

        !-------------------------------------------------------------------------------------
        ! calculate number of possible hoppings for each spin UP configuration
        ! Coeffs_Hubbard_nConnectedUP(i1) = number of configurations that can be reached by hoping
        ! Coeffs_Hubbard_hoppUP(i1,i) = configurations that can be reached from configuration i1 (index i)
        ! Coeffs_Hubbard_weightUP(i1,i) = matrix element for the hopping
        !-------------------------------------------------------------------------------------

        do iCoeff = 1, Coeffs_Hubbard_nCoeffsUP
          bits = Coeffs_Hubbard_bitcodesUP(iCoeff)

          tbits = ExciteBits([i1], [j1], bits)
          if (tbits .eq. 0) cycle

          Coeffs_Hubbard_nConnectedUP(iCoeff) = Coeffs_Hubbard_nConnectedUP(iCoeff) + 1
          count = Coeffs_Hubbard_nConnectedUP(iCoeff)

          ! save the configuration that can be reached
          Coeffs_Hubbard_hoppUP(count, iCoeff) = ConfigurationsFromIndexUP(abs(tbits))
          ! Use the precalculated one-body Hamiltonian matrix element
          Coeffs_Hubbard_weightUP(count, iCoeff) = real(Method_Mb_OrbBased_h1(i1, j1) * sign(1_I64, tbits), kind=R64)

        end do

        !-------------------------------------------------------------------------------------
        ! repeat for spin down
        !-------------------------------------------------------------------------------------

        do iCoeff = 1, Coeffs_Hubbard_nCoeffsDN
          bits = Coeffs_Hubbard_bitcodesDN(iCoeff)

          tbits = ExciteBits([i1], [j1], bits)
          if (tbits .eq. 0) cycle

          Coeffs_Hubbard_NConnectedDN(iCoeff) = Coeffs_Hubbard_NConnectedDN(iCoeff) + 1
          count = Coeffs_Hubbard_NConnectedDN(iCoeff)

          ! save the configuration that can be reached
          Coeffs_Hubbard_hoppDN(count, iCoeff) = ConfigurationsFromIndexDN(abs(tbits))
          ! Use the precalculated one-body Hamiltonian matrix element
          Coeffs_Hubbard_weightDN(count, iCoeff) = real(Method_Mb_OrbBased_h1(i1, j1) * sign(1_I64, tbits), kind=R64)

        end do
      end do
    end do

    do iCoeff = 1, Coeffs_Hubbard_nCoeffsUP

      call Combinatorics_sortintegerarray(Coeffs_Hubbard_hoppUP(:, iCoeff), Coeffs_Hubbard_nConnectedUP(iCoeff), permut_=permut)

      upTmp(:) = Coeffs_Hubbard_weightUP(:, iCoeff)

      do k = 1, Coeffs_Hubbard_nConnectedUP(iCoeff)
        Coeffs_Hubbard_weightUP(k, iCoeff) = upTmp(permut(k))
      end do

    end do

    do iCoeff = 1, Coeffs_Hubbard_nCoeffsDN

      call Combinatorics_sortintegerarray(Coeffs_Hubbard_hoppDN(:, iCoeff), Coeffs_Hubbard_NConnectedDN(iCoeff), permut_=permut)

      dnTmp(:) = Coeffs_Hubbard_weightDN(:, iCoeff)

      do k = 1, Coeffs_Hubbard_NConnectedDN(iCoeff)
        Coeffs_Hubbard_weightDN(k, iCoeff) = dnTmp(permut(k))
      end do

    end do

    !-------------------------------------------------------------------------------------
    ! Calculate the onsite interaction for each pair of up/down configurations
    ! by finding doubly occupied sites (bitwise AND of up and down configurations)
    !-------------------------------------------------------------------------------------

    do iUp = 1, Coeffs_Hubbard_nCoeffsUP
      do iDn = 1, Coeffs_Hubbard_nCoeffsDN
        ! Perform bitwise AND to find doubly occupied sites
        interactionBits = iand(Coeffs_Hubbard_bitcodesUP(iUp), Coeffs_Hubbard_bitcodesDN(iDn))

        ! Calculate interaction value based on the doubly occupied sites
        interactionStrength = 0.0_R64

        ! For each bit position (site), check if it's doubly occupied
        do iSite = 1, Grid_nPoints
          bitpos = iSite - 1

          ! If the bit is set (site is doubly occupied)
          if (btest(interactionBits, bitpos)) then
            ! Add the onsite interaction from h2 matrix for this site
            onSite = 2.0_R64 * real(Method_Mb_OrbBased_h2(iSite, iSite + Grid_nPoints, iSite, iSite + Grid_nPoints), kind=R64)
            interactionStrength = interactionStrength + onSite
          end if
        end do

        ! Store the interaction strength for this pair of configurations
        Coeffs_Hubbard_interactionValues(interactionBits) = interactionStrength

      end do

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function ExciteBits(creates, destroys, bits) result(res)

    integer(I64)              :: res
    integer(I32), intent(in), contiguous  :: creates(:)
    integer(I32), intent(in), contiguous  :: destroys(:)
    integer(I64), intent(in) :: bits
    integer(I32) :: sgn, iDestroy, iCreate, i, j

    res = bits
    sgn = 1

    do i = 1, size(destroys)
      iDestroy = destroys(i)

      if (.not. btest(res, iDestroy - 1)) then
        res = 0
        return
      end if

      do j = 1, iDestroy - 1
        if (btest(res, j - 1)) sgn = (-1) * sgn
      end do

      res = ibclr(res, iDestroy - 1)

    end do

    do i = 1, size(creates)
      iCreate = creates(size(creates) - i + 1)

      if (btest(res, iCreate - 1)) then
        res = 0
        return
      end if

      do j = 1, iCreate - 1
        if (btest(res, j - 1)) sgn = (-1) * sgn
      end do

      res = ibset(res, iCreate - 1)

    end do

    res = sgn * res

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyExcitation(coeffs, creates, destroys, bt)
    use M_Coeffs
    use M_Method_Mb

    complex(R64), intent(inout), contiguous  :: coeffs(:)
    integer(I32), intent(in), contiguous  :: creates(:)
    integer(I32), intent(in), contiguous  :: destroys(:)
    integer(I32), intent(in) :: bt

    complex(R64), allocatable  :: coeffsTmp(:)

    integer(I32) :: iCoeff, iCoeffNew, configurationBtNew
    integer(I64) :: bits, tbits
    real(R64)  :: factor
    integer(I32) :: configurations(Method_Mb_nBodyTypes)

    allocate (coeffsTmp, source=coeffs)

    coeffs(:) = 0.0_R64

    do iCoeff = 1, Coeffs_nCoeffs

      call Coeffs_ConfigurationsFromIndex(configurations, iCoeff)

      bits = Coeffs_Hubbard_bitcodesUP(configurations(bt))
      tbits = ExciteBits(creates, destroys, bits)
      if (tbits .eq. 0) cycle

      configurationBtNew = ConfigurationsFromIndexUP(abs(tbits))
      factor = sign(1.0_R64, real(tbits, kind=R64))

      configurations(bt) = configurationBtNew

      call Coeffs_IndexFromConfigurations(iCoeffNew, configurations)

      coeffs(iCoeffNew) = factor * coeffsTmp(iCoeff)

    end do

  end subroutine

end submodule
