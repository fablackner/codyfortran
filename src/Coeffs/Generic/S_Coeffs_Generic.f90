! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Implementation submodule for the Generic coefficient backend.
!>
!> This submodule contains the actual algorithms for:
!> - `ApplyH1FillRdm1`: One-body Hamiltonian application with optional 1-RDM
!> - `ApplyH2FillRdm2`: Two-body Hamiltonian application with optional 2-RDM
!> - `ApplyExcitation`: Ladder operator application (a†_i a_j sequences)
!> - `IndexFromConfigurations` / `ConfigurationsFromIndex`: Index mapping
!>
!> # Implementation Strategy
!>
!> **Hamiltonian Application:**
!> Uses sparse connectivity graphs from ConfigList. For each configuration:
!> 1. Loop over body types
!> 2. Use precomputed `singles` (for H1) or `doubles` (for H2 intra-body)
!> 3. For H2 inter-body: combine singles from different body types
!> 4. Accumulate contributions with OpenMP reduction
!>
!> **RDM Construction:**
!> Simultaneously builds RDM while applying H. The linearized `orbCode` packs
!> orbital indices into a single integer for cache-efficient accumulation.
!>
!> **Index Convention:**
!> A connection (orbCode(i,j), factor, jC) of configuration iC stores
!> factor·|jC⟩ = a†_i a_j |iC⟩, so the accumulation loops naturally act with
!> the transposed index order. The packing helpers (`FillH1Lin`, `FillH2Lin`,
!> `FillRdm1FromRdm1Lin`, `FillRdm2FromRdm2Lin`) translate between this
!> internal order and the public convention of M_Coeffs:
!>   h1(i,j) = ⟨i|h|j⟩,  rdm1(i,j) = ⟨a†_j a_i⟩,
!>   rdm2(i1,i2,j1,j2) = ⟨a†_j1 a†_j2 a_i2 a_i1⟩.
!> The transposition is invisible for real orbitals and Hamiltonians but
!> essential for complex ones (e.g., velocity-gauge laser coupling).
!>
!> **Total-M sector (coeffs.generic.totalM):**
!> With declared orbital m (method.mb.orbBased.mOfOrbs) every configuration
!> of a body type has a definite M, and h1, h2 conserve the total M. Setup
!> then lists the configurations of the declared total M (`sector`) and sorts
!> the connections of the ConfigList tables of every body type in place by
!> the change of M they cause (`sectorBt`: buckets of the singles per change,
!> the doubles that conserve M; no other reader depends on their order). The
!> Hamiltonian loops run over the sector and pair a single of change d on one
!> body type with a single of change -d on the other, so no connection leaves
!> the sector: exact for vectors in the sector, which Coeffs_ProjectOnSector
!> enforces once at initialization. Storage stays full length, and on it the
!> Hamiltonian acts as H_sector ⊕ 0 (see M_Coeffs). Without the
!> key the sector holds all configurations and one bucket all connections
!> (in their order of ConfigList), so the loops are those of the full space.
!>
!> **Excitation Application:**
!> Exploits tensor-product structure: only the affected body-type slice changes.
!> Uses a reusable buffer (`coeffsBuffer`) to avoid repeated allocation.
!>
!> # Data Structures
!>
!> - `h1Lin(nO²)`: Linearized 1-body Hamiltonian for fast lookup
!> - `h2Lin(nO⁴)`: Linearized 2-body Hamiltonian (symmetry-adapted storage)
!> - `rdm1Lin/rdm2Lin`: Thread-local accumulators (linearized)
!> - `coeffsBuffer`: Module-level buffer for excitation operations
!> - `sector`, `sectorBt`: CI sector and the buckets of its connections (see above)
submodule(M_Coeffs_Generic) S_Coeffs_Generic
  use M_ConfigList, only: T_ConfigList_Excitation

  implicit none

  !> Per body type: M of every configuration and the buckets of the sorted
  !> ConfigList tables (SortIntoMBuckets): the singles of configuration iC
  !> with change d of M are k = singlesFirst(d, iC) .. singlesFirst(d + 1, iC)
  !> - 1, d = -dMMax..dMMax; the doubles that conserve M are
  !> k = doublesFirst(0, iC) .. doublesFirst(1, iC) - 1.
  type :: T_SectorBt
    integer(I32), allocatable :: mOfConfig(:)
    integer(I32), allocatable :: singlesFirst(:, :)
    integer(I32), allocatable :: doublesFirst(:, :)
  end type

!=============================================================================
! module-scope persistent data
!=============================================================================

  !> Reusable buffer for ApplyExcitation to avoid repeated allocation
  complex(R64), allocatable, save :: coeffsBuffer(:)

  !> coeffs.generic.totalM declared, and its value (0 if not)
  logical, save      :: sectorQ
  integer(I32), save :: totalM

  !> Largest change of M of a single (0 without sector: one bucket)
  integer(I32), save :: dMMax

  !> Indices of the configurations of total M equal to totalM (all without sector)
  integer(I32), allocatable, save :: sector(:)

  !> M of the configurations and buckets of the connections per body type
  type(T_SectorBt), allocatable, save :: sectorBt(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Bind generic backend procedures and compute total CI dimension.
!>
!> Total coefficient count is the product over all body types:
!>   nCoeffs = Π_{bt} configList(bt)%nConfigurations
  module subroutine Coeffs_Generic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Method_Mb_OrbBased
    use M_Coeffs
    use M_ConfigList

    integer(I32) :: i

    call Say_Fabricate("coeffs.generic")

    !------------------------------------
    ! Compute total CI dimension
    !------------------------------------

    Coeffs_nCoeffs = 1
    do i = 1, size(configList)
      Coeffs_nCoeffs = Coeffs_nCoeffs * configList(i) % e % nConfigurations
    end do

    !------------------------------------
    ! Bind procedure pointers
    !------------------------------------

    Coeffs_ApplyH1FillRdm1 => ApplyH1FillRdm1
    Coeffs_ApplyH2FillRdm2 => ApplyH2FillRdm2
    Coeffs_ApplyExcitation => ApplyExcitation
    Coeffs_IndexFromConfigurations => IndexFromConfigurations
    Coeffs_ConfigurationsFromIndex => ConfigurationsFromIndex
    Coeffs_Setup => Setup

    !------------------------------------
    ! declared total-M sector
    !------------------------------------

    sectorQ = Json_GetExistence("coeffs.generic.totalM")
    totalM = 0
    if (sectorQ) then
      totalM = Json_Get("coeffs.generic.totalM", 0)
      if (.not. Method_Mb_OrbBased_definiteMQ) error stop "coeffs.generic.totalM requires method.mb.orbBased.mOfOrbs"
      Coeffs_ProjectOnSector => ProjectOnSector
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Build the sector and sort the ConfigList tables into its buckets. The M of
!> a configuration is the sum of the declared m over its occupations
!> n_i (a†_i a_i |C⟩ = n_i |C⟩); 0 without sector, so that the sector holds
!> all configurations and the buckets all connections.
  subroutine Setup
    use M_Utils_Say
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Coeffs
    use M_ConfigList

    type(T_SectorBt), allocatable :: bts(:)
    logical, allocatable :: inSectorQ(:)
    integer(I32) :: configurations(Method_Mb_nBodyTypes)
    integer(I32) :: bt, iC, iCNew, i, shift, iCoeff
    real(R64) :: occupation

    call Say_Setup("coeffs.generic")

    dMMax = 0
    if (sectorQ) dMMax = maxval(Method_Mb_OrbBased_mOfOrb) - minval(Method_Mb_OrbBased_mOfOrb)

    allocate (bts(Method_Mb_nBodyTypes))
    do bt = 1, Method_Mb_nBodyTypes
      associate (e => configList(bt) % e)
        allocate (bts(bt) % mOfConfig(e % nConfigurations), source=0)
        if (sectorQ) then
          shift = sum(Method_Mb_OrbBased_nOrbs(1:bt - 1))
          do iC = 1, e % nConfigurations
            do i = 1, Method_Mb_OrbBased_nOrbs(bt)
              call e % ExciteConfiguration(iCNew, occupation, [i], [i], iC)
              if (iCNew .ne. 0) bts(bt) % mOfConfig(iC) = bts(bt) % mOfConfig(iC) + &
                                                         nint(occupation) * Method_Mb_OrbBased_mOfOrb(shift + i)
            end do
          end do
        end if

        call SortIntoMBuckets(bts(bt) % singlesFirst, e % singles, bts(bt) % mOfConfig, dMMax)
        call SortIntoMBuckets(bts(bt) % doublesFirst, e % doubles, bts(bt) % mOfConfig, 0)
      end associate
    end do
    call move_alloc(bts, sectorBt)

    allocate (inSectorQ(Coeffs_nCoeffs))
    do iCoeff = 1, Coeffs_nCoeffs
      call ConfigurationsFromIndex(configurations, iCoeff)
      inSectorQ(iCoeff) = sum([(sectorBt(bt) % mOfConfig(configurations(bt)), bt=1, Method_Mb_nBodyTypes)]) .eq. totalM
    end do
    sector = pack([(iCoeff, iCoeff=1, Coeffs_nCoeffs)], inSectorQ)

    if (sectorQ) write (*, '(A,I0,A,I0,A,I0,A)') " coeffs.generic.totalM = ", totalM, ": sector of ", &
      size(sector), " of ", Coeffs_nCoeffs, " configurations"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Sort the connections of every configuration of `table` in place (stably)
!> by the change d = mOfConfig(target) - mOfConfig(iC) of M they cause, into
!> the buckets d = -dMax..dMax: first(d, iC) .. first(d + 1, iC) - 1.
!> Connections of larger change follow the last bucket and are not used; with
!> all d = 0 (no sector) the order stays unchanged.
  subroutine SortIntoMBuckets(first, table, mOfConfig, dMax)
    integer(I32), intent(out), allocatable       :: first(:, :)
    type(T_ConfigList_Excitation), intent(inout) :: table
    integer(I32), intent(in)                     :: mOfConfig(:)
    integer(I32), intent(in)                     :: dMax

    integer(I32), allocatable :: dOf(:), order(:)
    integer(I32) :: iC, n, d, k

    allocate (first(-dMax:dMax + 1, size(table % nConnected)))

    do iC = 1, size(table % nConnected)
      n = table % nConnected(iC)
      dOf = mOfConfig(table % excitedC(1:n, iC)) - mOfConfig(iC)
      where (abs(dOf) > dMax) dOf = dMax + 1
      order = [(pack([(k, k=1, n)], dOf .eq. d), d=-dMax, dMax + 1)]
      do d = -dMax, dMax + 1
        first(d, iC) = 1 + count(dOf < d)
      end do
      table % excitedC(1:n, iC) = table % excitedC(order, iC)
      table % factor(1:n, iC) = table % factor(order, iC)
      table % orbCode(1:n, iC) = table % orbCode(order, iC)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Zero the coefficients outside the sector; error stop if more than 1e-8 of
!> the norm would be removed.
  subroutine ProjectOnSector(coeffs)
    complex(R64), intent(inout), contiguous :: coeffs(:)

    real(R64), parameter :: tolerance = 1.0e-8_R64

    logical, allocatable :: outsideQ(:)
    real(R64) :: norm, normRemoved

    allocate (outsideQ(size(coeffs)))
    outsideQ(:) = .true.
    outsideQ(sector) = .false.

    norm = sum(abs(coeffs)**2)
    normRemoved = sum(abs(coeffs)**2, mask=outsideQ)
    if (normRemoved > tolerance * norm) then
      write (*, '(A,I0,A,ES10.3)') " fraction of the norm outside total M = ", totalM, ": ", normRemoved / norm
      error stop "Coeffs_ProjectOnSector: coefficients do not lie in the sector of coeffs.generic.totalM"
    end if

    where (outsideQ) coeffs = 0.0_R64

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Apply one-body Hamiltonian H1 and optionally compute 1-RDM.
!>
!> The one-body Hamiltonian is:
!>   Ĥ₁ = Σ_{pq} h¹_{pq} a†_p a_q
!>
!> This routine uses precomputed `singles` connectivity from ConfigList:
!> for each configuration iCoeff, iterate over allowed single excitations
!> factor·|jCoeff⟩ = a†_p a_q |iCoeff⟩ and accumulate:
!>   dCoeffs(iCoeff) += h1(q,p) × factor × coeffs(jCoeff)
!>
!> If `rdm1_` is requested, simultaneously accumulates:
!>   rdm1(q,p) += factor × conj(coeffs(jCoeff)) × coeffs(iCoeff)
!> (the transposition is done by the packing helpers, see Index Convention)
!> Rows and connections are those of the sector (singles of bucket 0).
!>
!> @param[in]    coeffs   CI coefficient vector
!> @param[out]   rdm1_    (optional) 1-RDM output, allocated if present
!> @param[inout] dCoeffs_ (optional) accumulator for H1|Ψ⟩
!> @param[in]    h1_      (optional) explicit 1-body matrix; uses internal if absent
  subroutine ApplyH1FillRdm1(coeffs, rdm1_, dCoeffs_, h1_)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Coeffs
    use M_ConfigList

    complex(R64), intent(in), contiguous, target                 :: coeffs(:)
    complex(R64), intent(out), allocatable, optional             :: rdm1_(:, :)
    complex(R64), intent(inout), contiguous, target, optional    :: dCoeffs_(:)
    complex(R64), intent(in), contiguous, optional               :: h1_(:, :)

    integer(I32) :: iBt1C, iS, iCoeff, jC
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

    !$omp parallel default(shared) private(iS, iCoeff, iTmp, iBt1, iBt1C, configurations, k1, orbCode, factor, jC, rdm1LinThread)

    if (present(rdm1_)) allocate (rdm1LinThread(nO**2))
    if (present(rdm1_)) rdm1LinThread(:) = 0.0_R64

    !$omp do
    do iS = 1, size(sector)

      iCoeff = sector(iS)
      iTmp = 1

      call Coeffs_ConfigurationsFromIndex(configurations, iCoeff)

      do iBt1 = 1, Method_Mb_nBodyTypes

        ! Extract configuration index for this body type
        iBt1C = configurations(iBt1)

        ! Loop over connected configurations via single excitations that conserve M
        do k1 = sectorBt(iBt1) % singlesFirst(0, iBt1C), sectorBt(iBt1) % singlesFirst(1, iBt1C) - 1

          orbCode = configList(iBt1) % e % singles % orbCode(k1, iBt1C)
          factor = configList(iBt1) % e % singles % factor(k1, iBt1C)
          jC = configList(iBt1) % e % singles % excitedC(k1, iBt1C)

          ! Transform body-type-local config change to global linear index
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
!> Apply two-body Hamiltonian H2 and optionally compute 2-RDM.
!>
!> The two-body Hamiltonian is:
!>   Ĥ₂ = Σ_{pqrs} h²_{pqrs} a†_p a†_q a_s a_r,  h²_{pqrs} = ½⟨pq|W|rs⟩
!> (index pairs are translated by the packing helpers, see Index Convention)
!>
!> Two interaction types are handled:
!> 1. **Intra-body** (same body type): uses `doubles` connectivity
!> 2. **Inter-body** (different body types): combines `singles` from both
!>
!> The inter-body term is non-trivial because excitations on different body
!> types commute, allowing factorization as (single on bt1) ⊗ (single on bt2).
!>
!> Rows and connections are those of the sector: doubles that conserve M, and
!> singles of change d on bt1 paired with singles of change -d on bt2.
!>
!> @param[in]    coeffs   CI coefficient vector
!> @param[out]   rdm2_    (optional) 2-RDM output, allocated if present
!> @param[inout] dCoeffs_ (optional) accumulator for H2|Ψ⟩
!> @param[in]    h2_      (optional) explicit 2-body tensor
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
    integer(I32) :: iS, iCoeff, jC, jCTmp
    integer(I32) :: k1, k2, d, iBt1, iBt2, iTmp, jTmp, orbCode, orbCodeTmp
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

    !$omp parallel default(shared) private(iS, iCoeff, iTmp, iBt1, iBt1C, configurations, k1, orbCode, factor, jC, orbCodeTmp, &
    !$omp                                  factorTmp, jCTmp, jTmp, iBt2, iBt2C, k2, d, rdm2LinThreads)

    if (present(rdm2_)) allocate (rdm2LinThreads(nO**4))
    if (present(rdm2_)) rdm2LinThreads(:) = 0.0_R64

    !$omp do
    do iS = 1, size(sector)

      iCoeff = sector(iS)
      iTmp = 1

      call Coeffs_ConfigurationsFromIndex(configurations, iCoeff)

      do iBt1 = 1, Method_Mb_nBodyTypes

        iBt1C = configurations(iBt1)

        !----------------------------------------------------------------------
        ! Intra-body interaction: double excitation within same body type
        !----------------------------------------------------------------------
        do k1 = sectorBt(iBt1) % doublesFirst(0, iBt1C), sectorBt(iBt1) % doublesFirst(1, iBt1C) - 1

          orbCode = configList(iBt1) % e % doubles % orbCode(k1, iBt1C)
          factor = configList(iBt1) % e % doubles % factor(k1, iBt1C)
          jC = configList(iBt1) % e % doubles % excitedC(k1, iBt1C)

          ! Map local config change to global index
          jC = (jC - iBt1C)
          jC = iCoeff + jC * iTmp

          if (present(dCoeffs_)) dCoeffs_(iCoeff) = dCoeffs_(iCoeff) + factor * h2Lin(orbCode) * coeffs(jC)
          if (present(rdm2_)) rdm2LinThreads(orbCode) = rdm2LinThreads(orbCode) + factor * conjg(coeffs(jC)) * coeffs(iCoeff)

        end do

        if (iBt1 .eq. Method_Mb_nBodyTypes) cycle

        !----------------------------------------------------------------------
        ! Inter-body interaction: combine singles from different body types,
        ! of change d of M on iBt1 and -d on iBt2
        !----------------------------------------------------------------------
        do d = -dMMax, dMMax
          do k1 = sectorBt(iBt1) % singlesFirst(d, iBt1C), sectorBt(iBt1) % singlesFirst(d + 1, iBt1C) - 1

            orbCodeTmp = configList(iBt1) % e % singles % orbCode(k1, iBt1C)
            orbCodeTmp = (orbCodeTmp - 1) * nO * nO

            factorTmp = configList(iBt1) % e % singles % factor(k1, iBt1C)
            jCTmp = configList(iBt1) % e % singles % excitedC(k1, iBt1C)

            jCTmp = (jCTmp - iBt1C)
            jCTmp = iCoeff + jCTmp * iTmp

            jTmp = iTmp * configList(iBt1) % e % nConfigurations

            do iBt2 = iBt1 + 1, Method_Mb_nBodyTypes

              iBt2C = configurations(iBt2)

              do k2 = sectorBt(iBt2) % singlesFirst(-d, iBt2C), sectorBt(iBt2) % singlesFirst(-d + 1, iBt2C) - 1

                orbCode = orbCodeTmp + configList(iBt2) % e % singles % orbCode(k2, iBt2C)
                factor = factorTmp * configList(iBt2) % e % singles % factor(k2, iBt2C)
                jC = configList(iBt2) % e % singles % excitedC(k2, iBt2C)

                jC = (jC - iBt2C)
                jC = jCTmp + jC * jTmp

                if (present(dCoeffs_)) dCoeffs_(iCoeff) = dCoeffs_(iCoeff) + factor * h2Lin(orbCode) * coeffs(jC)
                if (present(rdm2_)) rdm2LinThreads(orbCode) = rdm2LinThreads(orbCode) + factor * conjg(coeffs(jC)) * coeffs(iCoeff)

              end do

              jTmp = jTmp * configList(iBt2) % e % nConfigurations
            end do

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
!> Unpack linearized 1-RDM into standard 2D matrix form.
!>
!> The accumulator at `orbCode = (i-1)*nO + j` holds ⟨a†_i a_j⟩, which is
!> rdm1(j,i) in the public convention rdm1(i,j) = ⟨a†_j a_i⟩.
  subroutine FillRdm1FromRdm1Lin(rdm1, rdm1Lin)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: rdm1(:, :)
    complex(R64), intent(in), contiguous   :: rdm1Lin(:)

    integer(I32) :: i1, j1
    integer(I32) :: nO, orbCode

    nO = Method_Mb_OrbBased_nOrbsSum

    do j1 = 1, nO
      do i1 = 1, nO

        orbCode = 1
        orbCode = (orbCode - 1) * nO + i1
        orbCode = (orbCode - 1) * nO + j1

        rdm1(j1, i1) = rdm1Lin(orbCode)

      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Unpack linearized 2-RDM into standard 4D tensor form with symmetry expansion.
!>
!> The linearized form stores only unique elements (i1 ≤ i2, j1 ≤ j2).
!> This routine expands to all permutations with proper (anti)symmetry factors
!> based on particle statistics:
!> - Fermions: rdm2(i1,i2,j1,j2) = -rdm2(i2,i1,j1,j2) = -rdm2(i1,i2,j2,j1) = +rdm2(i2,i1,j2,j1)
!> - Bosons:   rdm2(i1,i2,j1,j2) = +rdm2(i2,i1,j1,j2) = +rdm2(i1,i2,j2,j1) = +rdm2(i2,i1,j2,j1)
!>
!> Inter-body terms (different body types) use antisymmetric convention via
!> Klein transformation for consistent handling of distinguishable particles.
!>
!> The accumulators hold ⟨a†_i1 a†_i2 a_j2 a_j1⟩; the expansion writes them
!> directly with creator and annihilator index pairs swapped, i.e. in the public
!> convention rdm2(i1,i2,j1,j2) = ⟨a†_j1 a†_j2 a_i2 a_i1⟩.
  subroutine FillRdm2FromRdm2Lin(rdm2, rdm2Lin)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm2Lin(:)

    integer(I32) :: i1, i2, j1, j2, i1bt, i2bt, j1bt, j2bt
    integer(I32) :: nO, orbCode

    nO = Method_Mb_OrbBased_nOrbsSum

    rdm2(:, :, :, :) = 0.0_R64

    do j2 = 1, nO
      do j1 = 1, j2
        do i2 = 1, nO
          do i1 = 1, i2

            j2bt = Method_Mb_OrbBased_bodyTypeOfOrb(j2)
            j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)
            i2bt = Method_Mb_OrbBased_bodyTypeOfOrb(i2)
            i1bt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)

            ! Skip invalid body-type combinations
            if (i1bt .eq. i2bt .and. j1bt .ne. j2bt) cycle
            if (i1bt .ne. i2bt .and. j1bt .eq. j2bt) cycle

            orbCode = 1
            orbCode = (orbCode - 1) * nO + i1
            orbCode = (orbCode - 1) * nO + j1
            orbCode = (orbCode - 1) * nO + i2
            orbCode = (orbCode - 1) * nO + j2

            if (i1bt .eq. i2bt) then
              ! Same body type: use actual statistics

              if (Method_Mb_bodyStatistics(i1bt) .eq. 'f') then
                ! Fermionic antisymmetry
                rdm2(j1, j2, i1, i2) = rdm2(j1, j2, i1, i2) + rdm2Lin(orbCode)
                rdm2(j1, j2, i2, i1) = rdm2(j1, j2, i2, i1) - rdm2Lin(orbCode)
                rdm2(j2, j1, i1, i2) = rdm2(j2, j1, i1, i2) - rdm2Lin(orbCode)
                rdm2(j2, j1, i2, i1) = rdm2(j2, j1, i2, i1) + rdm2Lin(orbCode)

              else if (Method_Mb_bodyStatistics(i1bt) .eq. 'b') then
                ! Bosonic symmetry (with diagonal correction)
                rdm2(j1, j2, i1, i2) = rdm2(j1, j2, i1, i2) + rdm2Lin(orbCode)
                rdm2(j1, j2, i2, i1) = rdm2(j1, j2, i2, i1) + rdm2Lin(orbCode)
                rdm2(j2, j1, i1, i2) = rdm2(j2, j1, i1, i2) + rdm2Lin(orbCode)
                rdm2(j2, j1, i2, i1) = rdm2(j2, j1, i2, i1) + rdm2Lin(orbCode)

                if (i1 .eq. i2) rdm2(j1, j2, i1, i2) = rdm2(j1, j2, i1, i2) / 2.0_R64
                if (j1 .eq. j2) rdm2(j1, j2, i1, i2) = rdm2(j1, j2, i1, i2) / 2.0_R64

              else
                error stop "Method_Mb_bodyStatistics not implemented"
              end if

            else
              ! Different body types: use antisymmetric convention (Klein transform)
              rdm2(j1, j2, i1, i2) = rdm2(j1, j2, i1, i2) + rdm2Lin(orbCode)
              rdm2(j1, j2, i2, i1) = rdm2(j1, j2, i2, i1) - rdm2Lin(orbCode)
              rdm2(j2, j1, i1, i2) = rdm2(j2, j1, i1, i2) - rdm2Lin(orbCode)
              rdm2(j2, j1, i2, i1) = rdm2(j2, j1, i2, i1) + rdm2Lin(orbCode)

            end if

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Pack 2D 1-body Hamiltonian into linearized form for fast lookup.
!>
!> Uses same encoding as `FillRdm1FromRdm1Lin`: orbCode = (i-1)*nO + j. The
!> connection a†_i a_j of configuration iC couples to ⟨iC|Ĥ₁|jC⟩ through
!> h1(j,i), so the transposed element is stored.
  subroutine FillH1Lin(h1Lin, h1)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: h1Lin(:)
    complex(R64), intent(in), contiguous   :: h1(:, :)

    integer(I32) :: i1, j1
    integer(I32) :: nO, orbCode

    nO = Method_Mb_OrbBased_nOrbsSum

    do j1 = 1, nO
      do i1 = 1, nO

        orbCode = 1
        orbCode = (orbCode - 1) * nO + i1
        orbCode = (orbCode - 1) * nO + j1

        h1Lin(orbCode) = h1(j1, i1)

      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Pack 4D 2-body Hamiltonian into linearized form with symmetry folding.
!>
!> Exploits permutation symmetry to store only unique elements (i1 ≤ i2, j1 ≤ j2).
!> The stored value incorporates the (anti)symmetry sum so that the lookup
!> during H2 application returns the correctly weighted matrix element.
!>
!> Same body type with fermionic statistics uses antisymmetrized form:
!>   h2Lin = h2(j1,j2,i1,i2) - h2(j2,j1,i1,i2) - h2(j1,j2,i2,i1) + h2(j2,j1,i2,i1)
!>
!> Same body type with bosonic statistics uses symmetrized form (with diagonal fix).
!> Different body types use antisymmetric convention via Klein transformation.
!>
!> As for `FillH1Lin`, the connections act with the transposed index order, so
!> h2 is read with its creator and annihilator index pairs swapped.
  subroutine FillH2Lin(h2Lin, h2)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: h2Lin(:)
    complex(R64), intent(in), contiguous   :: h2(:, :, :, :)

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

            ! Skip invalid body-type combinations
            if (i1bt .eq. i2bt .and. j1bt .ne. j2bt) cycle
            if (i1bt .ne. i2bt .and. j1bt .eq. j2bt) cycle

            orbCode = 1
            orbCode = (orbCode - 1) * nO + i1
            orbCode = (orbCode - 1) * nO + j1
            orbCode = (orbCode - 1) * nO + i2
            orbCode = (orbCode - 1) * nO + j2

            if (i1bt .eq. i2bt) then

              if (Method_Mb_bodyStatistics(i1bt) .eq. 'f') then
                ! Fermionic: antisymmetrize (vanishes on diagonal)
                h2Lin(orbCode) = (h2(j1, j2, i1, i2) - h2(j2, j1, i1, i2) - h2(j1, j2, i2, i1) + h2(j2, j1, i2, i1))
                if (i1 .eq. i2) h2Lin(orbCode) = 0
                if (j1 .eq. j2) h2Lin(orbCode) = 0

              else if (Method_Mb_bodyStatistics(i1bt) .eq. 'b') then
                ! Bosonic: symmetrize (with diagonal correction)
                h2Lin(orbCode) = (h2(j1, j2, i1, i2) + h2(j2, j1, i1, i2) + h2(j1, j2, i2, i1) + h2(j2, j1, i2, i1))

                if (i1 .eq. i2) h2Lin(orbCode) = h2Lin(orbCode) / 2
                if (j1 .eq. j2) h2Lin(orbCode) = h2Lin(orbCode) / 2

              else
                error stop "Method_Mb_bodyStatistics not implemented"
              end if

            else
              ! Different body types: antisymmetric convention
              h2Lin(orbCode) = (h2(j1, j2, i1, i2) - h2(j2, j1, i1, i2) - h2(j1, j2, i2, i1) + h2(j2, j1, i2, i1))
              if (i1 .eq. i2) h2Lin(orbCode) = 0
              if (j1 .eq. j2) h2Lin(orbCode) = 0

            end if

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Apply excitation operators a†_{creates} a_{destroys} for a specific body type.
!>
!> Transforms the CI vector by applying the operator sequence:
!>   a†_{creates(1)} a†_{creates(2)} ... a_{destroys(n)} ... a_{destroys(1)}
!>
!> Exploits tensor-product structure: only the bt-th factor of the direct product
!> is modified. For each configuration in body type `bt`, determine target config
!> and phase via `configList(bt)%ExciteConfiguration`, then scatter values using
!> precomputed strides.
!>
!> Uses a module-level buffer (`coeffsBuffer`) for efficiency; resized if needed.
!>
!> @param[inout] coeffs   CI coefficients (modified in-place)
!> @param[in]    creates  Orbital indices for creation operators
!> @param[in]    destroys Orbital indices for annihilation operators
!> @param[in]    bt       Body type index for this excitation
  subroutine ApplyExcitation(coeffs, creates, destroys, bt)
    use M_Coeffs
    use M_Method_Mb
    use M_ConfigList

    complex(R64), intent(inout), contiguous  :: coeffs(:)
    integer(I32), intent(in), contiguous  :: creates(:)
    integer(I32), intent(in), contiguous  :: destroys(:)
    integer(I32), intent(in) :: bt

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

    ! Compute stride: product of nConfigurations for body types < bt
    strideBt = 1
    do ibt_loop = 1, bt - 1
      strideBt = strideBt * configList(ibt_loop) % e % nConfigurations
    end do

    nConfigsBt = configList(bt) % e % nConfigurations
    stepSize = strideBt * nConfigsBt
    nUpper = Coeffs_nCoeffs / stepSize

    ! Iterate over configurations of the target body type
    do c = 1, nConfigsBt

      call configList(bt) % e % ExciteConfiguration(cNew, factor, creates, destroys, c)
      if (cNew .eq. 0) cycle

      offsetOld = (c - 1) * strideBt
      offsetNew = (cNew - 1) * strideBt

      ! Scatter to all tensor blocks (vectorized copy with factor)
      do k = 0, nUpper - 1
        loopIndex = k * stepSize
        coeffs(loopIndex + offsetNew + 1:loopIndex + offsetNew + strideBt) = &
          factor * coeffsBuffer(loopIndex + offsetOld + 1:loopIndex + offsetOld + strideBt)
      end do

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Convert linear CI index to configuration tuple (inverse of IndexFromConfigurations).
!>
!> Uses multi-base decomposition: iCoeff-1 = Σ_bt (C_bt-1) × stride_bt
!> where stride_bt = Π_{bt' < bt} nConfigurations(bt').
!>
!> @param[out] configurations  Configuration indices for each body type (1-based)
!> @param[in]  iCoeff          Linear CI index (1-based)
  pure subroutine ConfigurationsFromIndex(configurations, iCoeff)
    use M_ConfigList

    integer(I32), intent(out), contiguous :: configurations(:)
    integer(I32), intent(in) :: iCoeff

    integer(I32) :: ibt, iTmp

    iTmp = iCoeff - 1
    do ibt = 1, size(configurations)
      configurations(ibt) = mod(iTmp, configList(ibt) % e % nConfigurations) + 1
      iTmp = iTmp / configList(ibt) % e % nConfigurations
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Convert configuration tuple to linear CI index.
!>
!> Computes: iCoeff = 1 + Σ_bt (configurations(bt)-1) × stride_bt
!>
!> @param[out] iCoeff          Linear CI index (1-based)
!> @param[in]  configurations  Configuration indices for each body type (1-based)
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
