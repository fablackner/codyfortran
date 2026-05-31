! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_ConfigList_AllActive_Fermionic.f90
!> @brief Fermionic Fock-state enumeration and excitation logic.
!>
!> @details
!> ## Fermionic Configuration Encoding
!>
!> Configurations are represented as bit-patterns stored in `codeFromConfig`:
!> ```
!>   bit i = 1  ⟺  orbital i is occupied
!>   bit i = 0  ⟺  orbital i is empty (hole)
!> ```
!>
!> For N fermions in M orbitals, there are C(M,N) configurations. The bit-pattern
!> representation enables O(1) occupation queries via `btest(code, i-1)`.
!>
!> ## Anti-commutation and Sign Factor
!>
!> Fermionic creation/annihilation operators satisfy:
!> ```
!>   {a_i, a†_j} = δ_ij,   {a_i, a_j} = {a†_i, a†_j} = 0
!> ```
!>
!> When applying a_j to a configuration, we must count occupied orbitals to the
!> left of j to determine the phase:
!> ```
!>   a_j |n_1, ..., n_M⟩ = (-1)^{Σ_{k<j} n_k} × n_j × |..., n_j-1, ...⟩
!> ```
!>
!> The returned `factor` is ±1 encoding this parity.
!>
!> ## Pauli Exclusion
!>
!> - Creating on occupied orbital → iCNew = 0 (forbidden)
!> - Annihilating empty orbital  → iCNew = 0 (forbidden)
!>
!> ## Configuration Index Mapping
!>
!> Linear configuration index iC ∈ [1, C(M,N)] uses lexicographic ordering of
!> combinations without repetition. The utility `Combinatorics_indexOfCombiNoRepeat`
!> converts occupation list → index.
submodule(M_ConfigList_AllActive_Fermionic) S_ConfigList_AllActive_Fermionic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Allocate a fermionic all-active configuration element.
!>
!> @param[out] e     Polymorphic output (allocated to T_ConfigList_E_AllActive_Fermionic)
!> @param[in]  path  JSON path (unused here, read by Fabricate)
  module subroutine ConfigList_E_AllActive_Fermionic_Allocate(e, path)
    use M_Utils_UnusedVariables

    class(T_ConfigList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_ConfigList_E_AllActive_Fermionic :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Fabricate fermionic element: read params and compute nConfigurations.
!>
!> @details
!> Reads from JSON:
!> - `bodyTarget`: which particle species (default 1)
!> - `nExcitations`: max excitation rank (default = nBodies = full CI)
!>
!> Computes nConfigurations as sum over excitation ranks:
!> ```
!>   nConfigs = Σ_{k=0}^{nE} C(M-N, k) × C(N, k)
!> ```
!> where M = nOrbitals, N = nBodies, and we're summing over k particle-hole
!> excitations from the reference.
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_SfGslLib
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_ConfigList

    class(T_ConfigList_E_AllActive_Fermionic), intent(inout) :: this

    integer(I32) :: bt, nOBt, nBBt, nE, i

    call Say_Fabricate(this % path//".fermionic")

    !------------------------------------
    ! read JSON parameters
    !------------------------------------

    this % bodyTarget = Json_Get("bodyTarget", 1, path_=this % path//".fermionic")
    bt = this % bodyTarget

    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)

    this % nExcitations = Json_Get("nExcitations", nBBt, path_=this % path//".fermionic")
    nE = this % nExcitations

    !------------------------------------
    ! compute number of configurations
    !------------------------------------
    ! Sum over excitation ranks: each k-fold excitation contributes
    ! C(nOrbs - nBodies, k) × C(nBodies, k) configurations

    this % nConfigurations = 0
    do i = 0, nE
      this % nConfigurations = this % nConfigurations + SfGslLib_Binomial(nOBt - nBBt, i) * SfGslLib_Binomial(nBBt, i)
    end do

    !------------------------------------
    ! validate statistics
    !------------------------------------

    if (Method_Mb_bodyStatistics(this % bodyTarget) .ne. 'f') error stop "bodyTarget not fermionic"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Setup: build codeFromConfig bit-pattern mapping.
!>
!> @details
!> ## Algorithm
!>
!> 1. Generate all C(M,N) combinations of N orbitals from M using
!>    `CombinationGslLib_CombiNoRepeat` → `i1(1:N, 1:nConfigs)`
!>
!> 2. Convert each combination to bit-pattern:
!>    ```
!>    codeFromConfig(iC) = Σ_j 2^{i1(j,iC) - 1}
!>    ```
!>
!> 3. Verify round-trip: `indexOfCombiNoRepeat(i1(:,iC)) == iC` for all configs.
!>
!> ## Thread Safety
!>
!> After Setup completes, codeFromConfig is read-only.
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_Combinatorics
    use M_Utils_CombinationGslLib
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive_Fermionic), intent(inout) :: this

    integer(I32) :: iC, j, bt, nConfigurations, nBBt, nOBt
    integer(I32), allocatable :: i1(:, :)

    call Say_Setup(this % path//".fermionic")

    allocate (this % codeFromConfig(this % nConfigurations))

    bt = this % bodyTarget
    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)
    nConfigurations = this % nConfigurations

    allocate (i1(nBBt, nConfigurations))

    !---------------------------------------------------------------------------
    ! i1(i, iC) is the orbital index of the i-th particle in config iC
    ! Generated in lexicographic order by GSL combination iterator
    !---------------------------------------------------------------------------

    call CombinationGslLib_CombiNoRepeat(i1, nOBt)

    !---------------------------------------------------------------------------
    ! Convert occupation lists to bit-patterns
    ! bit k is set ⟺ orbital k+1 is occupied
    !---------------------------------------------------------------------------

    this % codeFromConfig = 0
    do j = 1, nBBt
      do iC = 1, nConfigurations
        this % codeFromConfig(iC) = ibset(this % codeFromConfig(iC), i1(j, iC) - 1)
      end do
    end do

    !---------------------------------------------------------------------------
    ! Verify bijection between config index and occupation list
    !---------------------------------------------------------------------------

    do iC = 1, nConfigurations
      if ((Combinatorics_indexOfCombiNoRepeat(nOBt, i1(:, iC)) - iC) .ne. 0) error stop "index calc failed"
    end do

    deallocate (i1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Apply fermionic excitation: a†_{creates} a_{destroys} |iC⟩
!>
!> @param[in]  this     Configuration element
!> @param[out] iCNew    Resulting configuration index (0 if forbidden)
!> @param[out] factor   Sign factor ±1 from anti-commutation
!> @param[in]  creates  Orbital indices for creation operators
!> @param[in]  destroys Orbital indices for annihilation operators
!> @param[in]  iC       Input configuration index
!>
!> @details
!> ## Algorithm
!>
!> 1. Decode bit-pattern → occupation array
!>
!> 2. Apply annihilations (left to right):
!>    - Check orbital occupied (else return iCNew=0)
!>    - Accumulate sign: (-1)^{count of occupied orbitals to the left}
!>    - Decrement occupation
!>
!> 3. Apply creations (right to left, reversed order):
!>    - Check orbital empty (else return iCNew=0, Pauli exclusion)
!>    - Accumulate sign
!>    - Increment occupation
!>
!> 4. Convert final occupation → config index via `indexOfCombiNoRepeat`
!>
!> ## Why Reversed Creation Order?
!>
!> Physics convention: a†_1 a†_2 applies a†_2 first, then a†_1.
!> We reverse the creates array to match this operator ordering.
  module subroutine ExciteConfiguration(this, iCNew, factor, creates, destroys, iC)
    use M_Utils_Combinatorics
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive_Fermionic), intent(in) :: this
    integer(I32), intent(out)        :: iCNew
    real(R64), intent(out)           :: factor
    integer(I32), intent(in), contiguous          :: creates(:)
    integer(I32), intent(in), contiguous          :: destroys(:)
    integer(I32), intent(in)         :: iC

    integer(I32) :: nBBt, nOBt
    integer(I32) :: iDestroy, iCreate, i, j, sum

    integer(I64) :: code

    integer(I32), allocatable :: i1(:)
    integer(I32), allocatable :: occupation(:)

    nOBt = Method_Mb_OrbBased_nOrbs(this % bodyTarget)
    nBBt = Method_Mb_nBodies(this % bodyTarget)

    factor = 1.0_R64

    allocate (i1(nBBt))
    allocate (occupation(nOBt))

    !---------------------------------------------------------------------------
    ! Decode bit-pattern to occupation array
    !---------------------------------------------------------------------------

    code = this % codeFromConfig(iC)

    occupation = 0
    do i = 1, nOBt
      if (btest(code, i - 1)) then
        occupation(i) = 1
      end if
    end do

    !---------------------------------------------------------------------------
    ! Apply annihilation operators (left to right)
    !---------------------------------------------------------------------------

    do i = 1, size(destroys)
      iDestroy = destroys(i)

      ! Check orbital is occupied
      if (occupation(iDestroy) .eq. 0) then
        iCNew = 0
        return
      end if

      ! Accumulate fermionic sign: count occupied orbitals to the left
      do j = 1, iDestroy - 1
        if (occupation(j) .eq. 1) factor = (-1) * factor
      end do

      occupation(iDestroy) = occupation(iDestroy) - 1

    end do

    !---------------------------------------------------------------------------
    ! Apply creation operators (reversed order: rightmost first)
    !---------------------------------------------------------------------------

    do i = 1, size(creates)
      iCreate = creates(size(creates) - i + 1)

      ! Pauli exclusion: cannot create on occupied orbital
      if (occupation(iCreate) .eq. 1) then
        iCNew = 0
        return
      end if

      ! Accumulate fermionic sign
      do j = 1, iCreate - 1
        if (occupation(j) .eq. 1) factor = (-1) * factor
      end do

      occupation(iCreate) = occupation(iCreate) + 1

    end do

    !---------------------------------------------------------------------------
    ! Convert occupation array back to orbital list and config index
    !---------------------------------------------------------------------------

    i = 0
    sum = 1
    do while (sum <= nBBt)

      i = i + 1
      if (occupation(i) .eq. 0) cycle

      i1(sum) = i
      sum = sum + occupation(i)
    end do

    iCNew = Combinatorics_indexOfCombiNoRepeat(nOBt, i1(:))

  end subroutine

end submodule
