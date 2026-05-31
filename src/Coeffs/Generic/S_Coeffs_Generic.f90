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
submodule(M_Coeffs_Generic) S_Coeffs_Generic

  implicit none

!=============================================================================
! module-scope persistent data
!=============================================================================

  !> Reusable buffer for ApplyExcitation to avoid repeated allocation
  complex(R64), allocatable, save :: coeffsBuffer(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Bind generic backend procedures and compute total CI dimension.
!>
!> Total coefficient count is the product over all body types:
!>   nCoeffs = Π_{bt} configList(bt)%nConfigurations
  module subroutine Coeffs_Generic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
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

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Apply one-body Hamiltonian H1 and optionally compute 1-RDM.
!>
!> The one-body Hamiltonian is:
!>   Ĥ₁ = Σ_{pq} h¹_{pq} a†_p a_q
!>
!> This routine uses precomputed `singles` connectivity from ConfigList:
!> for each configuration iCoeff, iterate over allowed single excitations
!> (p→q) and accumulate:
!>   dCoeffs(iCoeff) += h1(p,q) × factor × coeffs(jCoeff)
!>
!> If `rdm1_` is requested, simultaneously accumulates:
!>   rdm1(p,q) += factor × conj(coeffs(jCoeff)) × coeffs(iCoeff)
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

        ! Extract configuration index for this body type
        iBt1C = configurations(iBt1)

        ! Loop over connected configurations via single excitations
        do k1 = 1, configList(iBt1) % e % singles % nConnected(iBt1C)

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
!>   Ĥ₂ = ½ Σ_{pqrs} h²_{pqrs} a†_p a†_q a_s a_r
!>
!> Two interaction types are handled:
!> 1. **Intra-body** (same body type): uses `doubles` connectivity
!> 2. **Inter-body** (different body types): combines `singles` from both
!>
!> The inter-body term is non-trivial because excitations on different body
!> types commute, allowing factorization as (single on bt1) ⊗ (single on bt2).
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

        iBt1C = configurations(iBt1)

        !----------------------------------------------------------------------
        ! Intra-body interaction: double excitation within same body type
        !----------------------------------------------------------------------
        do k1 = 1, configList(iBt1) % e % doubles % nConnected(iBt1C)

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
        ! Inter-body interaction: combine singles from different body types
        !----------------------------------------------------------------------
        do k1 = 1, configList(iBt1) % e % singles % nConnected(iBt1C)

          orbCodeTmp = configList(iBt1) % e % singles % orbCode(k1, iBt1C)
          orbCodeTmp = (orbCodeTmp - 1) * nO * nO

          factorTmp = configList(iBt1) % e % singles % factor(k1, iBt1C)
          jCTmp = configList(iBt1) % e % singles % excitedC(k1, iBt1C)

          jCTmp = (jCTmp - iBt1C)
          jCTmp = iCoeff + jCTmp * iTmp

          jTmp = iTmp * configList(iBt1) % e % nConfigurations

          do iBt2 = iBt1 + 1, Method_Mb_nBodyTypes

            iBt2C = configurations(iBt2)

            do k2 = 1, configList(iBt2) % e % singles % nConnected(iBt2C)

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
!> Maps `orbCode = (i-1)*nO + j` back to rdm1(i,j).
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
                rdm2_(i1, i2, j1, j2) = rdm2_(i1, i2, j1, j2) + rdm2Lin(orbCode)
                rdm2_(i2, i1, j1, j2) = rdm2_(i2, i1, j1, j2) - rdm2Lin(orbCode)
                rdm2_(i1, i2, j2, j1) = rdm2_(i1, i2, j2, j1) - rdm2Lin(orbCode)
                rdm2_(i2, i1, j2, j1) = rdm2_(i2, i1, j2, j1) + rdm2Lin(orbCode)

              else if (Method_Mb_bodyStatistics(i1bt) .eq. 'b') then
                ! Bosonic symmetry (with diagonal correction)
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
              ! Different body types: use antisymmetric convention (Klein transform)
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
!> Pack 2D 1-body Hamiltonian into linearized form for fast lookup.
!>
!> Uses same encoding as `FillRdm1FromRdm1Lin`: orbCode = (i-1)*nO + j.
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
!> Pack 4D 2-body Hamiltonian into linearized form with symmetry folding.
!>
!> Exploits permutation symmetry to store only unique elements (i1 ≤ i2, j1 ≤ j2).
!> The stored value incorporates the (anti)symmetry sum so that the lookup
!> during H2 application returns the correctly weighted matrix element.
!>
!> Same body type with fermionic statistics uses antisymmetrized form:
!>   h2Lin = h2(i1,i2,j1,j2) - h2(i1,i2,j2,j1) - h2(i2,i1,j1,j2) + h2(i2,i1,j2,j1)
!>
!> Same body type with bosonic statistics uses symmetrized form (with diagonal fix).
!> Different body types use antisymmetric convention via Klein transformation.
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
                h2Lin(orbCode) = (h2_(i1, i2, j1, j2) - h2_(i1, i2, j2, j1) - h2_(i2, i1, j1, j2) + h2_(i2, i1, j2, j1))
                if (i1 .eq. i2) h2Lin(orbCode) = 0
                if (j1 .eq. j2) h2Lin(orbCode) = 0

              else if (Method_Mb_bodyStatistics(i1bt) .eq. 'b') then
                ! Bosonic: symmetrize (with diagonal correction)
                h2Lin(orbCode) = (h2_(i1, i2, j1, j2) + h2_(i1, i2, j2, j1) + h2_(i2, i1, j1, j2) + h2_(i2, i1, j2, j1))

                if (i1 .eq. i2) h2Lin(orbCode) = h2Lin(orbCode) / 2
                if (j1 .eq. j2) h2Lin(orbCode) = h2Lin(orbCode) / 2

              else
                error stop "Method_Mb_bodyStatistics not implemented"
              end if

            else
              ! Different body types: antisymmetric convention
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
