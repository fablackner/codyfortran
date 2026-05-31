! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Implementation submodule for the Hubbard coefficient backend.
!>
!> This submodule contains:
!> - `Coeffs_Hubbard_Fabricate`: Backend selection and pointer binding
!> - `Setup`: Build bit-encoded configurations and hopping graphs
!> - `SetupConfigurations`: Enumerate configs via GSL combination generator
!> - `SetupHamiltonian`: Precompute hopping connectivity and interaction values
!> - `ExciteBits`: Apply creation/annihilation to bit pattern with sign tracking
!> - `ApplyExcitation`: Transform CI vector via ladder operators
!>
!> # Implementation Details
!>
!> **Configuration enumeration:** Uses GSL's combination iterator to generate
!> all C(n,k) ways to place k particles in n sites. Each combination is
!> converted to a bit pattern for O(1) operations.
!>
!> **Hopping graph construction:** For each (source, target) site pair with
!> nonzero h1 matrix element, check if the hop is valid (source occupied,
!> target empty). If so, store target config index and weight (including
!> fermionic sign from bit counting).
!>
!> **Fermionic sign:** When moving particle from site i to site j, count
!> occupied sites between i and j. Sign = (-1)^count. Implemented via
!> sequential bit testing in `ExciteBits`.
submodule(M_Coeffs_Hubbard) S_Coeffs_Hubbard

  implicit none

  !=============================================================================
  ! local data: inverse lookup tables for bit patterns
  !=============================================================================

  !> Inverse lookup: bit pattern → config index (spin up)
  integer(I32), allocatable :: ConfigurationsFromIndexUP(:)

  !> Inverse lookup: bit pattern → config index (spin down)
  integer(I32), allocatable :: ConfigurationsFromIndexDN(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Bind Hubbard backend procedures and select spin-symmetry variant.
!>
!> Computes per-spin configuration counts via binomial coefficients:
!>   nCoeffsUP = C(nOrbs, nBodiesUp)
!>   nCoeffsDN = C(nOrbs, nBodiesDn)
!>
!> Then delegates to one of the spin-symmetry sub-backends based on JSON.
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
    ! Compute per-spin configuration counts
    !------------------------------------

    Coeffs_Hubbard_nCoeffsUP = SfGslLib_Binomial(Method_Mb_OrbBased_nOrbs(1), Method_Mb_nBodies(1))
    Coeffs_Hubbard_nCoeffsDN = SfGslLib_Binomial(Method_Mb_OrbBased_nOrbs(2), Method_Mb_nBodies(2))

    !------------------------------------
    ! Common procedure pointers
    !------------------------------------

    Coeffs_Setup => Setup
    Coeffs_ApplyExcitation => ApplyExcitation

    !------------------------------------
    ! Select spin-symmetry variant from JSON
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
!> Initialize Hubbard backend data structures.
!>
!> Must be called after Fabricate. Builds:
!> 1. Bit-encoded configuration arrays
!> 2. Inverse lookup tables (bit pattern → index)
!> 3. Hopping connectivity graphs
!> 4. Precomputed interaction values
  subroutine Setup
    use M_Utils_Say

    call Say_Setup("coeffs.hubbard")

    call SetupConfigurations
    call SetupHamiltonian

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Enumerate all configurations and build bit patterns + inverse lookup.
!>
!> Uses GSL combination generator to enumerate all ways to place N particles
!> in M sites. Each combination is converted to a bit pattern where bit j=1
!> means site j is occupied.
!>
!> Inverse lookup tables allow O(1) conversion: bit pattern → config index.
  subroutine SetupConfigurations
    use M_Utils_CombinationGslLib
    use M_Grid
    use M_Method_Mb

    integer(I32) :: iCoeff, j
    integer(I32), allocatable :: collect(:, :)

    allocate (Coeffs_Hubbard_bitcodesUP(Coeffs_Hubbard_nCoeffsUP))
    allocate (Coeffs_Hubbard_bitcodesDN(Coeffs_Hubbard_nCoeffsDN))

    ! Inverse lookup: bit pattern → config index (needs 2^nSites entries)
    allocate (ConfigurationsFromIndexUP(0:2**(Grid_nPoints)))
    allocate (ConfigurationsFromIndexDN(0:2**(Grid_nPoints)))

    !-------------------------------------------------------------------------------------
    ! Spin-up configurations
    !-------------------------------------------------------------------------------------

    allocate (collect(Method_Mb_nBodies(1), Coeffs_Hubbard_nCoeffsUP))

    ! collect(i, iCoeff) = position of i-th particle in configuration iCoeff
    call CombinationGslLib_CombiNoRepeat(collect, Grid_nPoints)

    ! Convert particle positions to bit patterns
    Coeffs_Hubbard_bitcodesUP(:) = 0
    do j = 1, Method_Mb_nBodies(1)
      do iCoeff = 1, Coeffs_Hubbard_nCoeffsUP
        Coeffs_Hubbard_bitcodesUP(iCoeff) = ibset(Coeffs_Hubbard_bitcodesUP(iCoeff), collect(j, iCoeff) - 1)
      end do
    end do

    ! Build inverse lookup table
    ConfigurationsFromIndexUP(:) = 0
    do iCoeff = 1, Coeffs_Hubbard_nCoeffsUP
      ConfigurationsFromIndexUP(Coeffs_Hubbard_bitcodesUP(iCoeff)) = iCoeff
    end do

    deallocate (collect)

    !-------------------------------------------------------------------------------------
    ! Spin-down configurations (same procedure)
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
!> Build hopping connectivity graphs and precompute interaction values.
!>
!> For each source configuration i and target site pair (i1,j1) with h1(i1,j1)≠0:
!> - Check if hop is valid: site i1 occupied and site j1 empty
!> - If valid, compute target config via `ExciteBits` (handles fermionic sign)
!> - Store in sparse format: hoppUP(k,i), weightUP(k,i), nConnectedUP(i)
!>
!> Interaction values are precomputed for all (up,dn) config pairs:
!>   interactionValues(bits) = U × (number of doubly occupied sites)
!> where bits = bitcodesUP(i) AND bitcodesDN(j).
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

    ! Validate prerequisites
    if (.not. allocated(Method_Mb_OrbBased_h1)) then
      error stop "Method_Mb_OrbBased_h1 is not allocated, please call Method_Mb_OrbBased_FillH1 first"
    end if

    if (.not. allocated(Method_Mb_OrbBased_h2)) then
      error stop "Method_Mb_OrbBased_h2 is not allocated, please call Method_Mb_OrbBased_FillH2 first"
    end if

    allocate (Coeffs_Hubbard_nConnectedUP(Coeffs_Hubbard_nCoeffsUP))
    allocate (Coeffs_Hubbard_NConnectedDN(Coeffs_Hubbard_nCoeffsDN))

    ! Count active dimensions for array sizing
    nDim = 0
    if (Grid_Lattice_xSize > 1) nDim = nDim + 1
    if (Grid_Lattice_ySize > 1) nDim = nDim + 1
    if (Grid_Lattice_zSize > 1) nDim = nDim + 1

    ! Max hops per config: nParticles × nDimensions × 2 (forward/backward)
    allocate (Coeffs_Hubbard_hoppUP(Method_Mb_nBodies(1) * nDim * 2, Coeffs_Hubbard_nCoeffsUP))
    allocate (Coeffs_Hubbard_hoppDN(Method_Mb_nBodies(2) * nDim * 2, Coeffs_Hubbard_nCoeffsDN))

    allocate (Coeffs_Hubbard_weightUP(Method_Mb_nBodies(1) * nDim * 2, Coeffs_Hubbard_nCoeffsUP))
    allocate (Coeffs_Hubbard_weightDN(Method_Mb_nBodies(2) * nDim * 2, Coeffs_Hubbard_nCoeffsDN))

    allocate (Coeffs_Hubbard_interactionValues(2**Grid_nPoints))

    allocate (upTmp(Method_Mb_nBodies(1) * nDim * 2))
    allocate (dnTmp(Method_Mb_nBodies(1) * nDim * 2))

    !-------------------------------------------------------------------------------------
    ! Build hopping connectivity for each (source, target) pair with nonzero h1
    !-------------------------------------------------------------------------------------

    Coeffs_Hubbard_nConnectedUP(:) = 0
    Coeffs_Hubbard_nConnectedDN(:) = 0

    Do i1 = 1, Grid_nPoints
      Do j1 = 1, Grid_nPoints
        if (abs(Method_Mb_OrbBased_h1(i1, j1)) < 1e-15_R64) cycle

        !---------------------------------------------------------------------------
        ! Spin-up hopping
        !---------------------------------------------------------------------------
        do iCoeff = 1, Coeffs_Hubbard_nCoeffsUP
          bits = Coeffs_Hubbard_bitcodesUP(iCoeff)

          tbits = ExciteBits([i1], [j1], bits)
          if (tbits .eq. 0) cycle

          Coeffs_Hubbard_nConnectedUP(iCoeff) = Coeffs_Hubbard_nConnectedUP(iCoeff) + 1
          count = Coeffs_Hubbard_nConnectedUP(iCoeff)

          Coeffs_Hubbard_hoppUP(count, iCoeff) = ConfigurationsFromIndexUP(abs(tbits))
          Coeffs_Hubbard_weightUP(count, iCoeff) = real(Method_Mb_OrbBased_h1(i1, j1) * sign(1_I64, tbits), kind=R64)
        end do

        !---------------------------------------------------------------------------
        ! Spin-down hopping (same procedure)
        !---------------------------------------------------------------------------
        do iCoeff = 1, Coeffs_Hubbard_nCoeffsDN
          bits = Coeffs_Hubbard_bitcodesDN(iCoeff)

          tbits = ExciteBits([i1], [j1], bits)
          if (tbits .eq. 0) cycle

          Coeffs_Hubbard_NConnectedDN(iCoeff) = Coeffs_Hubbard_NConnectedDN(iCoeff) + 1
          count = Coeffs_Hubbard_NConnectedDN(iCoeff)

          Coeffs_Hubbard_hoppDN(count, iCoeff) = ConfigurationsFromIndexDN(abs(tbits))
          Coeffs_Hubbard_weightDN(count, iCoeff) = real(Method_Mb_OrbBased_h1(i1, j1) * sign(1_I64, tbits), kind=R64)
        end do
      end do
    end do

    ! Sort hopping arrays for cache-friendly access during Hamiltonian application
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
    ! Precompute on-site interaction for all (up,dn) configuration pairs
    !-------------------------------------------------------------------------------------

    do iUp = 1, Coeffs_Hubbard_nCoeffsUP
      do iDn = 1, Coeffs_Hubbard_nCoeffsDN
        ! Doubly occupied sites = bitwise AND of up and down patterns
        interactionBits = iand(Coeffs_Hubbard_bitcodesUP(iUp), Coeffs_Hubbard_bitcodesDN(iDn))

        interactionStrength = 0.0_R64
        do iSite = 1, Grid_nPoints
          bitpos = iSite - 1
          if (btest(interactionBits, bitpos)) then
            onSite = 2.0_R64 * real(Method_Mb_OrbBased_h2(iSite, iSite + Grid_nPoints, iSite, iSite + Grid_nPoints), kind=R64)
            interactionStrength = interactionStrength + onSite
          end if
        end do

        Coeffs_Hubbard_interactionValues(interactionBits) = interactionStrength
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Apply ladder operators to a bit pattern, computing fermionic sign.
!>
!> Transforms `bits` by:
!> 1. Annihilating particles at sites in `destroys` (right to left order)
!> 2. Creating particles at sites in `creates` (left to right order)
!>
!> Returns 0 if any annihilation targets an empty site or creation targets
!> an occupied site (Pauli exclusion violation).
!>
!> The result encodes the fermionic sign in its sign: abs(result) is the
!> new bit pattern, sign(result) is ±1.
!>
!> **Sign calculation:** For each operator, count occupied sites to the left
!> (i.e., with smaller index). Each occupied site contributes a factor of -1.
!>
!> @param[in]  creates   Sites to create particles on (1-based)
!> @param[in]  destroys  Sites to annihilate particles from (1-based)
!> @param[in]  bits      Input bit pattern
!> @return     Signed result: |res| = new pattern, sign(res) = phase
  function ExciteBits(creates, destroys, bits) result(res)

    integer(I64)              :: res
    integer(I32), intent(in), contiguous  :: creates(:)
    integer(I32), intent(in), contiguous  :: destroys(:)
    integer(I64), intent(in) :: bits
    integer(I32) :: sgn, iDestroy, iCreate, i, j

    res = bits
    sgn = 1

    ! Apply annihilation operators (right to left)
    do i = 1, size(destroys)
      iDestroy = destroys(i)

      ! Check if site is occupied
      if (.not. btest(res, iDestroy - 1)) then
        res = 0
        return
      end if

      ! Count occupied sites to the left for fermionic sign
      do j = 1, iDestroy - 1
        if (btest(res, j - 1)) sgn = (-1) * sgn
      end do

      res = ibclr(res, iDestroy - 1)
    end do

    ! Apply creation operators (left to right, reverse order)
    do i = 1, size(creates)
      iCreate = creates(size(creates) - i + 1)

      ! Check if site is empty (Pauli exclusion)
      if (btest(res, iCreate - 1)) then
        res = 0
        return
      end if

      ! Count occupied sites to the left for fermionic sign
      do j = 1, iCreate - 1
        if (btest(res, j - 1)) sgn = (-1) * sgn
      end do

      res = ibset(res, iCreate - 1)
    end do

    res = sgn * res

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Apply ladder operators to CI vector for a specific body type.
!>
!> Transforms all CI coefficients by applying the operator sequence:
!>   a†_{creates(1)} a†_{creates(2)} ... a_{destroys(n)} ... a_{destroys(1)}
!> to the `bt`-th spin sector.
!>
!> Uses `ExciteBits` to compute target configuration and phase for each
!> source configuration, then scatters the result.
!>
!> @note Currently uses the spin-up inverse lookup table for all body types,
!>       which assumes equal numbers of up and down configurations. This is
!>       valid for half-filling but may need generalization for other cases.
!>
!> @param[inout] coeffs   CI coefficients (modified in-place)
!> @param[in]    creates  Orbital indices for creation operators
!> @param[in]    destroys Orbital indices for annihilation operators
!> @param[in]    bt       Body type index (1=up, 2=down)
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
