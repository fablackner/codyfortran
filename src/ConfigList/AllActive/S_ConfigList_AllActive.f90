! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_ConfigList_AllActive.f90
!> @brief Implementation of all-active configuration fabrication and connectivity setup.
!>
!> @details
!> This submodule provides:
!> 1. **Allocator dispatch** – routes to fermionic or bosonic implementations
!> 2. **Common setup orchestration** – calls level-2 setup then builds connectivity
!> 3. **Connectivity graph construction** – populates `singles` and `doubles`
!>
!> ## Connectivity Graph Algorithm
!>
!> For each configuration iC, we iterate over all possible excitations and call
!> `ExciteConfiguration`. Valid excitations (iCNew ≠ 0) are stored in sparse
!> column-major arrays. This is O(nConfigs × nOrbs^2) for singles and
!> O(nConfigs × nOrbs^4) for doubles, but only done once during Setup.
!>
!> ## OpenMP Parallelization
!>
!> The outer loop over configurations is parallelized. Each thread processes
!> independent configurations; no locking needed since writes are to disjoint
!> columns.
!>
!> ## orbCode Encoding
!>
!> Orbital indices are packed into a single integer for cache-friendly storage:
!> ```
!>   Singles:  orbCode = 1 * nO + i  then  orbCode * nO + j
!>             → decodes as a†_i a_j
!>
!>   Doubles:  orbCode = (((1 * nO + i1) * nO + j1) * nO + i2) * nO + j2
!>             → decodes as a†_i1 a†_i2 a_j2 a_j1
!> ```
!> The leading "1" ensures non-zero orbCode even for orbital index 0.
submodule(M_ConfigList_AllActive) S_ConfigList_AllActive

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Allocate concrete all-active element based on statistics in JSON.
!>
!> @param[out] e    Allocated configuration element (polymorphic output)
!> @param[in]  path JSON path to this element (e.g., "configList.allActive1")
!>
!> @details
!> Examines `path.fermionic` and `path.bosonic` to determine which concrete
!> type to allocate. Exactly one must be present.
  module subroutine ConfigList_AllActive_Allocate(e, path)
    use M_Utils_Json
    use M_Utils_Say
    use M_ConfigList_AllActive_Fermionic, only: ConfigList_E_AllActive_Fermionic_Allocate
    use M_ConfigList_AllActive_Bosonic, only: ConfigList_E_AllActive_Bosonic_Allocate

    class(T_ConfigList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    !------------------------------------
    ! branch on statistics
    !------------------------------------

    if (Json_GetExistence(path//".fermionic")) then
      call ConfigList_E_AllActive_Fermionic_Allocate(e, path//".fermionic")

    else if (Json_GetExistence(path//".bosonic")) then
      call ConfigList_E_AllActive_Bosonic_Allocate(e, path//".bosonic")

    else
      error stop path//". is missing one of: fermionic, bosonic"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Fabrication entry point for all-active elements.
!>
!> @details
!> Delegates immediately to statistics-specific FabricateLevel2. The common
!> all-active fabrication logic is minimal; most work happens at level-2.
  module subroutine Fabricate(this)
    use M_Utils_Say

    class(T_ConfigList_E_AllActive), intent(inout) :: this

    call Say_Fabricate(this % path)

    call this % FabricateLevel2()

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Setup entry point: builds codeFromConfig and connectivity graphs.
!>
!> @details
!> Execution order:
!> 1. `SetupLevel2` – statistics-specific: allocate codeFromConfig, build mapping
!> 2. `SetupSinglesData` – enumerate all single excitations, store connectivity
!> 3. `SetupDoublesData` – enumerate all double excitations, store connectivity
!>
!> After this, the element is fully initialized and thread-safe for read access.
  module subroutine Setup(this)
    use M_Utils_Say

    class(T_ConfigList_E_AllActive), intent(inout) :: this

    call Say_Setup(this % path)

    call this % SetupLevel2()

    call SetupSinglesData(this)
    call SetupDoublesData(this)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Build single-excitation connectivity graph.
!>
!> @details
!> ## Algorithm
!>
!> **Pass 1:** Count connections per configuration (OpenMP parallel).
!> ```
!>   for each iC:
!>     for j = 1..nOrbs, i = 1..nOrbs:
!>       if ExciteConfiguration(a†_i a_j, iC) succeeds:
!>         nConnected(iC) += 1
!> ```
!>
!> **Allocate:** Size arrays to maxConnections across all configs.
!>
!> **Pass 2:** Store (excitedC, factor, orbCode) for each valid excitation.
!>
!> ## Memory Layout
!>
!> Column-major: `excitedC(k, iC)` is the k-th target from config iC.
!> This gives cache-friendly iteration when processing one config at a time.
  subroutine SetupSinglesData(this)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive), intent(inout) :: this

    integer(I32) :: i1, j1, iC, iCNew
    integer(I32) :: nConnectedSinglesMax
    real(R64) :: factor
    integer(I32) :: k, bt, nOBt, nConfigurations, orbCode, shift, nO, nBBt

    bt = this % bodyTarget
    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)
    nO = Method_Mb_OrbBased_nOrbsSum
    shift = sum(Method_Mb_OrbBased_nOrbs(1:bt - 1))
    nConfigurations = this % nConfigurations

    allocate (this % singles % nConnected(nConfigurations))

    !---------------------------------------------------------------------------
    ! Pass 1: count valid connections per configuration
    !---------------------------------------------------------------------------
    !$omp parallel default(shared) private (iC, k, j1, i1, iCNew, factor)

    !$omp do
    do iC = 1, nConfigurations

      k = 0

      do j1 = 1, nOBt
        do i1 = 1, nOBt

          call this % ExciteConfiguration(iCNew, factor, [i1], [j1], iC)
          if (iCNew .eq. 0) cycle

          k = k + 1
        end do
      end do

      this % singles % nConnected(iC) = k

    end do
    !$omp end do

    !$omp end parallel

    nConnectedSinglesMax = maxval(this % singles % nConnected(:))

    allocate (this % singles % excitedC(nConnectedSinglesMax, nConfigurations))
    allocate (this % singles % orbCode(nConnectedSinglesMax, nConfigurations))
    allocate (this % singles % factor(nConnectedSinglesMax, nConfigurations))

    !---------------------------------------------------------------------------
    ! Pass 2: store connectivity data
    !---------------------------------------------------------------------------
    !$omp parallel default(shared) private (iC, k, j1, i1, iCNew, factor, orbCode)

    !$omp do
    do iC = 1, nConfigurations

      k = 0

      do j1 = 1, nOBt
        do i1 = 1, nOBt

          call this % ExciteConfiguration(iCNew, factor, [i1], [j1], iC)
          if (iCNew .eq. 0) cycle

          k = k + 1

          ! Encode orbital indices: orbCode = (i1 + shift) * nO + (j1 + shift)
          ! with leading 1 to avoid zero for i1=j1=0
          orbCode = 1
          orbCode = (orbCode - 1) * nO + (i1 + shift)
          orbCode = (orbCode - 1) * nO + (j1 + shift)

          this % singles % orbCode(k, iC) = orbCode
          this % singles % excitedC(k, iC) = iCNew
          this % singles % factor(k, iC) = factor

        end do
      end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Build double-excitation connectivity graph.
!>
!> @details
!> Same two-pass algorithm as SetupSinglesData, but iterates over all pairs
!> (i1 ≤ i2, j1 ≤ j2) to avoid double-counting symmetric excitations.
!>
!> ## Iteration Order
!>
!> ```
!>   for j2 = 1..nOrbs:
!>     for j1 = 1..j2:        ← triangular: j1 ≤ j2
!>       for i2 = 1..nOrbs:
!>         for i1 = 1..i2:    ← triangular: i1 ≤ i2
!>           apply a†_i1 a†_i2 a_j2 a_j1
!> ```
!>
!> This covers all unique double excitations for both fermions (where i1≠i2
!> by Pauli) and bosons (where i1=i2 is valid).
  subroutine SetupDoublesData(this)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive), intent(inout) :: this

    integer(I32) :: i1, i2, j1, j2, iC, iCNew
    integer(I32) :: k, bt, nOBt, nConfigurations, orbCode, shift, nO, nBBt
    integer(I32) :: nConnectedDoublesMax
    real(R64) :: factor

    bt = this % bodyTarget
    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)
    nO = Method_Mb_OrbBased_nOrbsSum
    shift = sum(Method_Mb_OrbBased_nOrbs(1:bt - 1))
    nConfigurations = this % nConfigurations

    allocate (this % doubles % nConnected(nConfigurations))

    !---------------------------------------------------------------------------
    ! Pass 1: count valid connections per configuration
    !---------------------------------------------------------------------------
    !$omp parallel default(shared) private (iC, k, j2, i2, j1, i1, iCNew, factor)

    !$omp do
    do iC = 1, nConfigurations

      k = 0

      do j2 = 1, nOBt
        do j1 = 1, j2

          do i2 = 1, nOBt
            do i1 = 1, i2

              call this % ExciteConfiguration(iCNew, factor, [i1, i2], [j1, j2], iC)
              if (iCNew .eq. 0) cycle

              k = k + 1
            end do
          end do
        end do
      end do

      this % doubles % nConnected(iC) = k

    end do
    !$omp end do

    !$omp end parallel

    nConnectedDoublesMax = maxval(this % doubles % nConnected(:))

    allocate (this % doubles % excitedC(nConnectedDoublesMax, nConfigurations))
    allocate (this % doubles % orbCode(nConnectedDoublesMax, nConfigurations))
    allocate (this % doubles % factor(nConnectedDoublesMax, nConfigurations))

    !---------------------------------------------------------------------------
    ! Pass 2: store connectivity data
    !---------------------------------------------------------------------------
    !$omp parallel default(shared) private (iC, k, j2, i2, j1, i1, iCNew, factor, orbCode)

    !$omp do
    do iC = 1, nConfigurations

      k = 0

      do j2 = 1, nOBt
        do j1 = 1, j2

          do i2 = 1, nOBt
            do i1 = 1, i2

              call this % ExciteConfiguration(iCNew, factor, [i1, i2], [j1, j2], iC)
              if (iCNew .eq. 0) cycle

              k = k + 1

              ! Encode 4 orbital indices into single integer
              orbCode = 1
              orbCode = (orbCode - 1) * nO + (i1 + shift)
              orbCode = (orbCode - 1) * nO + (j1 + shift)
              orbCode = (orbCode - 1) * nO + (i2 + shift)
              orbCode = (orbCode - 1) * nO + (j2 + shift)

              this % doubles % orbCode(k, iC) = orbCode
              this % doubles % excitedC(k, iC) = iCNew
              this % doubles % factor(k, iC) = factor

            end do
          end do
        end do
      end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine

end submodule
