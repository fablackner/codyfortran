! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_ConfigList_AllActive_Bosons.f90
!> @brief Bosonic Fock-state enumeration and excitation logic.
!>
!> @details
!> ## Bosonic Configuration Encoding
!>
!> Configurations are represented in base-(N+1) encoding stored in `codeFromConfig`:
!> ```
!>   digit i of code (in base N+1) = occupation of orbital i
!> ```
!>
!> For N bosons in M orbitals, there are C(M+N-1, N) configurations (stars-and-bars).
!> The encoding is normalized by dividing by N (since total occupation = N always).
!>
!> ## Ladder Factors
!>
!> Bosonic creation/annihilation operators satisfy:
!> ```
!>   [a_i, a†_j] = δ_ij,   [a_i, a_j] = [a†_i, a†_j] = 0
!> ```
!>
!> Unlike fermions, there's no sign factor. Instead, ladder prefactors arise:
!> ```
!>   a_i |n_i⟩ = √n_i |n_i - 1⟩
!>   a†_i |n_i⟩ = √(n_i + 1) |n_i + 1⟩
!> ```
!>
!> The returned `factor = √(product of all ladder prefactors)`.
!>
!> ## Configuration Index Mapping
!>
!> Linear configuration index iC ∈ [1, C(M+N-1,N)] uses lexicographic ordering of
!> combinations with repetition (multisets). The utility `Combinatorics_indexOfCombiWithRepeat`
!> converts occupation list → index.
!>
!> ## No Pauli Exclusion
!>
!> Multiple bosons can occupy the same orbital. Annihilating an empty orbital
!> (n_i = 0) returns iCNew = 0 (factor becomes 0).
submodule(M_ConfigList_AllActive_Bosonic) S_ConfigList_AllActive_Bosonic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Allocate a bosonic all-active configuration element.
!>
!> @param[out] e     Polymorphic output (allocated to T_ConfigList_E_AllActive_Bosonic)
!> @param[in]  path  JSON path (unused here, read by Fabricate)
  module subroutine ConfigList_E_AllActive_Bosonic_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_ConfigList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_ConfigList_E_AllActive_Bosonic :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Fabricate bosonic element: read params and compute nConfigurations.
!>
!> @details
!> Reads from JSON:
!> - `bodyTarget`: which particle species (default 1)
!> - `nExcitations`: max excitation rank (default = nBodies)
!>
!> Computes nConfigurations using stars-and-bars formula:
!> ```
!>   nConfigs = C(M + nE - 1, nE)
!> ```
!> where M = nOrbitals, nE = nExcitations.
!>
!> @note For bosons, nExcitations typically equals nBodies (full active space).
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_SfGslLib
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive_Bosonic), intent(inout) :: this

    integer(I32) :: bt, nOBt, nBBt, nE

    call Say_Fabricate(this % path//".bosonic")

    !------------------------------------
    ! read JSON parameters
    !------------------------------------

    this % bodyTarget = Json_Get("bodyTarget", 1, path_=this % path//".bosonic")
    bt = this % bodyTarget

    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)

    this % nExcitations = Json_Get("nExcitations", nBBt, path_=this % path//".bosonic")
    nE = this % nExcitations

    !------------------------------------
    ! compute number of configurations
    !------------------------------------
    ! Stars-and-bars: distribute N bosons among M orbitals

    this % nConfigurations = SfGslLib_Binomial(nOBt + nE - 1, nE)

    !------------------------------------
    ! validate statistics
    !------------------------------------

    if (Method_Mb_bodyStatistics(this % bodyTarget) .ne. 'b') error stop "bodyTarget not bosonic"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Setup: build codeFromConfig base-(N+1) encoding.
!>
!> @details
!> ## Algorithm
!>
!> 1. Generate all C(M+N-1,N) multisets of N orbitals from M using
!>    `MultisetGslLib_CombiWithRepeat` → `i1(1:N, 1:nConfigs)`
!>
!> 2. Convert each multiset to base-(N+1) encoding:
!>    ```
!>    codeFromConfig(iC) = Σ_j (N+1)^{i1(j,iC) - 1}
!>    ```
!>    This counts each orbital's contribution to the occupation.
!>
!> 3. Normalize by dividing by N (total occupation is always N).
!>
!> 4. Verify round-trip: `indexOfCombiWithRepeat(i1(:,iC)) == iC` for all configs.
!>
!> ## Why Base-(N+1)?
!>
!> With N bosons, any single orbital can have occupation 0..N. Base-(N+1)
!> uniquely encodes all possibilities. Dividing by N compresses the representation
!> since the sum of digits always equals N.
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_Combinatorics
    use M_Utils_MultisetGslLib
    use M_Utils_SfGslLib
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive_Bosonic), intent(inout) :: this

    integer(I32) :: iC, j, bt, nConfigurations, nOBt
    integer(I32), allocatable :: i1(:, :)

    integer(I32) :: nBBt

    call Say_Setup(this % path//".bosonic")

    allocate (this % codeFromConfig(this % nConfigurations))

    bt = this % bodyTarget
    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)
    nConfigurations = this % nConfigurations

    allocate (i1(nBBt, nConfigurations))

    !---------------------------------------------------------------------------
    ! i1(i, iC) is the orbital index of the i-th particle in config iC
    ! Generated in lexicographic order by GSL multiset iterator
    ! (with repetition allowed, so i1(1,iC) ≤ i1(2,iC) ≤ ... for sorted configs)
    !---------------------------------------------------------------------------

    call MultisetGslLib_CombiWithRepeat(i1, nOBt)

    !---------------------------------------------------------------------------
    ! Convert occupation lists to base-(N+1) encoding
    ! Each particle at orbital k contributes (N+1)^{k-1} to the code
    !---------------------------------------------------------------------------

    this % codeFromConfig(:) = 0
    do iC = 1, nConfigurations

      do j = 1, nBBt
        this % codeFromConfig(iC) = this % codeFromConfig(iC) + (nBBt + 1)**(i1(j, iC) - 1)
      end do

    end do

    !---------------------------------------------------------------------------
    ! Normalize by dividing by N
    ! Since total occupation = N, the raw code is divisible by N
    ! This compression reduces storage requirements
    !---------------------------------------------------------------------------

    this % codeFromConfig(:) = this % codeFromConfig(:) / nBBt

    !---------------------------------------------------------------------------
    ! Verify bijection between config index and occupation list
    !---------------------------------------------------------------------------

    do iC = 1, nConfigurations
      if ((Combinatorics_indexOfCombiWithRepeat(nOBt, i1(:, iC)) - iC) .ne. 0) error stop "index calc failed"
    end do

    deallocate (i1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Apply bosonic excitation: a†_{creates} a_{destroys} |iC⟩
!>
!> @param[in]  this     Configuration element
!> @param[out] iCNew    Resulting configuration index (0 if forbidden)
!> @param[out] factor   Ladder factor √(product of occupation factors)
!> @param[in]  creates  Orbital indices for creation operators
!> @param[in]  destroys Orbital indices for annihilation operators
!> @param[in]  iC       Input configuration index
!>
!> @details
!> ## Algorithm
!>
!> 1. Decode base-(N+1) code → occupation array (undo normalization first)
!>
!> 2. Apply annihilations:
!>    - Accumulate ladder factor: iFactor *= n_orb
!>    - If n_orb = 0, return iCNew = 0 (cannot annihilate vacuum)
!>    - Decrement occupation
!>
!> 3. Apply creations (reversed order to match operator convention):
!>    - Accumulate ladder factor: iFactor *= (n_orb + 1)
!>    - Increment occupation
!>
!> 4. Compute factor = √iFactor
!>
!> 5. Convert final occupation → config index via `indexOfCombiWithRepeat`
!>
!> ## Ladder Factor Accumulation
!>
!> For a_orb: multiply by current occupation (gives 0 if empty)
!> For a†_orb: multiply by (occupation + 1) before incrementing
!>
!> The final factor is the square root since we're working with
!> normalized states: a|n⟩ = √n |n-1⟩, a†|n⟩ = √(n+1) |n+1⟩.
  module subroutine ExciteConfiguration(this, iCNew, factor, creates, destroys, iC)
    use M_Utils_Combinatorics
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive_Bosonic), intent(in) :: this
    integer(I32), intent(out)        :: iCNew
    real(R64), intent(out)           :: factor
    integer(I32), intent(in), contiguous          :: creates(:)
    integer(I32), intent(in), contiguous          :: destroys(:)
    integer(I32), intent(in)         :: iC

    integer(I32) :: nBBt, nOBt
    integer(I32) :: iFactor, base
    integer(I32) :: orb, i
    integer(I32) :: sum

    integer(I64) :: code
    integer(I64) :: codeTmp

    integer(I32), allocatable :: i1(:)
    integer(I32), allocatable :: occupation(:)

    nOBt = Method_Mb_OrbBased_nOrbs(this % bodyTarget)
    nBBt = Method_Mb_nBodies(this % bodyTarget)

    base = nBBt + 1
    iFactor = 1

    allocate (i1(nBBt))
    allocate (occupation(nOBt))

    !---------------------------------------------------------------------------
    ! Decode base-(N+1) code to occupation array
    ! Undo normalization by multiplying by N first
    !---------------------------------------------------------------------------

    code = this % codeFromConfig(iC)
    codeTmp = code * nBBt

    do i = 1, nOBt
      occupation(i) = int(Mod(codeTmp, int(base, kind=I64)), kind=I32)
      codeTmp = codeTmp / base
    end do

    !---------------------------------------------------------------------------
    ! Apply annihilation operators
    !---------------------------------------------------------------------------

    do i = 1, size(destroys)
      orb = destroys(i)

      ! Accumulate ladder factor (becomes 0 if orbital empty)
      iFactor = iFactor * occupation(orb)

      if (iFactor .eq. 0) then
        iCNew = 0
        return
      end if

      occupation(orb) = occupation(orb) - 1

    end do

    !---------------------------------------------------------------------------
    ! Apply creation operators (reversed order for operator convention)
    !---------------------------------------------------------------------------

    do i = 1, size(creates)
      orb = creates(size(creates) - i + 1)

      ! Accumulate ladder factor: √(n+1) for creation
      iFactor = iFactor * (1 + occupation(orb))

      occupation(orb) = occupation(orb) + 1

    end do

    ! Take square root of accumulated integer factor
    factor = sqrt(real(iFactor, R64))

    !---------------------------------------------------------------------------
    ! Convert occupation array back to orbital list and config index
    !---------------------------------------------------------------------------

    i = 0
    sum = 1
    do while (sum <= nBBt)

      i = i + 1
      if (occupation(i) .eq. 0) cycle

      i1(sum:sum + occupation(i) - 1) = i
      sum = sum + occupation(i)
    end do

    iCNew = Combinatorics_indexOfCombiWithRepeat(nOBt, i1)

  end subroutine

end submodule
