! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Implementation submodule for M_Coeffs.
!>
!> Contains:
!> - `Coeffs_Fabricate`: JSON-driven backend selection (generic vs. Hubbard)
!> - Shared utility routines used by all backends:
!>   - `Normalize`:           L² normalization using BLAS
!>   - `ProjectOnSubspace`:   Gram-Schmidt projection onto orthogonal complement
!>   - `SaveCoeffs`:          Binary checkpoint output
!>   - `SaveTwoRdm`:          Compute and persist 2-RDM per body-type block
!>   - `FillRdm{1,2,3}Bt`:    Generic RDM construction via excitation application
!>   - `ApplyGroupedExcitations`: Helper to batch excitations by body type
!>
!> The body-type-resolved RDM routines (`FillRdm1Bt`, `FillRdm2Bt`, `FillRdm3Bt`)
!> use a brute-force loop over all orbital indices, applying ladder operators
!> via `Coeffs_ApplyExcitation` and computing expectation values. These are
!> primarily useful for diagnostics and validation; optimized backends may
!> provide direct implementations.
submodule(M_Coeffs) S_Coeffs

  implicit none

!=============================================================================
! local procedures
!=============================================================================

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Bind procedure pointers to the appropriate backend based on JSON config.
!>
!> Reads the `"coeffs"` section from global JSON and dispatches to either
!> `Coeffs_Generic_Fabricate` or `Coeffs_Hubbard_Fabricate`. Common routines
!> (normalization, projection, RDM fills, I/O) are assigned unconditionally
!> as they are backend-independent.
!>
!> **JSON keys:**
!> - `coeffs.generic` → M_Coeffs_Generic (tensor-product via ConfigList)
!> - `coeffs.hubbard` → M_Coeffs_Hubbard (bit-encoded lattice model)
  module subroutine Coeffs_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Coeffs_Generic
    use M_Coeffs_Hubbard

    call Say_Fabricate("coeffs")

    !------------------------------------
    ! Assign backend-independent shared procedures
    !------------------------------------

    Coeffs_Normalize => Normalize
    Coeffs_ProjectOnSubspace => ProjectOnSubspace
    Coeffs_SaveCoeffs => SaveCoeffs
    Coeffs_SaveTwoRdm => SaveTwoRdm
    Coeffs_FillRdm3Bt => FillRdm3Bt
    Coeffs_FillRdm2Bt => FillRdm2Bt
    Coeffs_FillRdm1Bt => FillRdm1Bt

    !------------------------------------
    ! Select backend from JSON
    !------------------------------------

    if (Json_GetExistence("coeffs.generic")) then
      call Coeffs_Generic_Fabricate

    else if (Json_GetExistence("coeffs.hubbard")) then
      call Coeffs_Hubbard_Fabricate

    else
      error stop "coeffs is missing one of: generic, hubbard"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Normalize CI vector to unit L² norm.
!>
!> Uses BLAS `dznrm2` (via BlasLib wrapper) for numerical stability on large
!> vectors. Modifies `coeffs` in-place.
  subroutine Normalize(coeffs)
    use M_Utils_BlasLib

    complex(R64), intent(inout), contiguous  :: coeffs(:)

    real(R64) :: norm

    norm = BlasLib_CalcNorm(coeffs)
    coeffs(:) = coeffs(:) / norm

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Gram-Schmidt projection: orthogonalize `dCoeffs` against `coeffs`.
!>
!> Computes:
!>   dCoeffs ← dCoeffs - ⟨coeffs|dCoeffs⟩ · coeffs
!>
!> This removes the component of `dCoeffs` parallel to `coeffs`, leaving only
!> the orthogonal complement. Assumes `coeffs` is normalized (or caller
!> accepts scaling).
  subroutine ProjectOnSubspace(dCoeffs, coeffs)
    complex(R64), intent(inout), contiguous  :: dCoeffs(:)
    complex(R64), intent(in), contiguous     :: coeffs(:)

    dCoeffs(:) = dCoeffs(:) - dot_product(coeffs, dCoeffs) * coeffs(:)

  end subroutine

!------------------------------------------------------------------------------
!> Group excitations by body type and apply sequentially.
!>
!> When computing RDM elements involving multiple body types (e.g., ⟨a†_α a†_β a_β a_α⟩),
!> excitations must be applied per body type. This helper:
!> 1. Finds unique body types in the excitation list
!> 2. Groups creation/annihilation indices by body type
!> 3. Calls `Coeffs_ApplyExcitation` once per group
!>
!> @param[inout] coeffsTmp  Working copy of CI vector (modified in-place)
!> @param[in]    iExc       Creation orbital indices (1..nExc)
!> @param[in]    jExc       Annihilation orbital indices (1..nExc)
!> @param[in]    btExc      Body type for each excitation (1..nExc)
!> @param[in]    nExc       Number of excitations
!------------------------------------------------------------------------------
  subroutine ApplyGroupedExcitations(coeffsTmp, iExc, jExc, btExc, nExc)
    complex(R64), intent(inout), contiguous :: coeffsTmp(:)
    integer(I32), intent(in) :: iExc(:)
    integer(I32), intent(in) :: jExc(:)
    integer(I32), intent(in) :: btExc(:)
    integer, intent(in) :: nExc

    integer(I32) :: uniqBts(nExc)
    integer(I32) :: iBuf(nExc), jBuf(nExc)
    integer      :: nGroups, g, k, cnt

    ! Collect unique body types
    nGroups = 0
    do k = 1, nExc
      if (.not. any(uniqBts(1:nGroups) .eq. btExc(k))) then
        nGroups = nGroups + 1
        uniqBts(nGroups) = btExc(k)
      end if
    end do

    ! Apply excitations grouped by body type
    do g = 1, nGroups
      cnt = 0
      do k = 1, nExc
        if (btExc(k) .eq. uniqBts(g)) then
          cnt = cnt + 1
          iBuf(cnt) = iExc(k)
          jBuf(cnt) = jExc(k)
        end if
      end do
      call Coeffs_ApplyExcitation(coeffsTmp, iBuf(1:cnt), jBuf(1:cnt), uniqBts(g))
    end do
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Save CI coefficients to binary file for checkpoint/restart.
!>
!> Writes to `coeffs.in` using the DataStorage module's raw binary format.
!> File can be loaded via CoeffsInit with `"load": {}` configuration.
  subroutine SaveCoeffs(coeffs)
    use M_Utils_DataStorage

    complex(R64), intent(in), contiguous  :: coeffs(:)

    call SaveData(coeffs, 'coeffs.in', storage_size(coeffs))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Compute full 2-RDM and save per body-type block.
!>
!> Produces files `rdm2BtXX_YY.in` for each body-type pair (XX, YY), where
!> the indices are zero-padded to 2 digits. These files can be loaded by
!> TwoRdmInit or used for external analysis.
!>
!> The 2-RDM is extracted from the full combined orbital space into blocks:
!>   rdm2Bt(i1, i2, j1, j2) with i1,j1 ∈ orbitals(bt1) and i2,j2 ∈ orbitals(bt2)
  subroutine SaveTwoRdm(coeffs)
    use M_Utils_DataStorage
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(in), contiguous  :: coeffs(:)

    complex(R64), allocatable :: rdm2(:, :, :, :)
    complex(R64), allocatable :: rdm2Bt(:, :, :, :)
    character(len=256) :: filename
    integer(I32)      :: ibt1, ibt2

    call Coeffs_ApplyH2FillRdm2(coeffs, rdm2_=rdm2)

    do ibt1 = 1, Method_Mb_nBodyTypes
      do ibt2 = 1, Method_Mb_nBodyTypes

        write (filename, '("rdm2Bt",I2.2,"_",I2.2,".in")') ibt1, ibt2

        rdm2Bt = rdm2(Method_Mb_OrbBased_nOrbsStart(ibt1):Method_Mb_OrbBased_nOrbsEnd(ibt1), &
                      Method_Mb_OrbBased_nOrbsStart(ibt2):Method_Mb_OrbBased_nOrbsEnd(ibt2), &
                      Method_Mb_OrbBased_nOrbsStart(ibt1):Method_Mb_OrbBased_nOrbsEnd(ibt1), &
                      Method_Mb_OrbBased_nOrbsStart(ibt2):Method_Mb_OrbBased_nOrbsEnd(ibt2))

        call SaveData(rdm2Bt, trim(filename), storage_size(rdm2Bt))

      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Compute 3-body RDM for specified body-type triplet.
!>
!> Computes:
!>   ρ(i1,i2,i3,j1,j2,j3) = ⟨Ψ| a†_{i1,bt1} a†_{i2,bt2} a†_{i3,bt3}
!>                                a_{j3,bt3} a_{j2,bt2} a_{j1,bt1} |Ψ⟩
!>
!> where orbital indices i,j run over the orbitals of the respective body type.
!>
!> **Complexity:** O(n₁² n₂² n₃² · nCoeffs) — use for small systems only.
!>
!> @param[out] rdm3Bt  6D array (n1×n2×n3×n1×n2×n3) allocated on output
!> @param[in]  coeffs  CI coefficient vector
!> @param[in]  bt1,bt2,bt3  Body type indices (1..nBodyTypes)
  subroutine FillRdm3Bt(rdm3Bt, coeffs, bt1, bt2, bt3)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm3Bt(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous, target     :: coeffs(:)
    integer(I32), intent(in) :: bt1
    integer(I32), intent(in) :: bt2
    integer(I32), intent(in) :: bt3

    complex(R64), allocatable  :: coeffsTmp(:)

    integer(I32) :: i1, i2, i3, j1, j2, j3

    if (.not. allocated(rdm3Bt)) allocate (rdm3Bt(Method_Mb_OrbBased_nOrbs(bt1), &
                                                  Method_Mb_OrbBased_nOrbs(bt2), &
                                                  Method_Mb_OrbBased_nOrbs(bt3), &
                                                  Method_Mb_OrbBased_nOrbs(bt1), &
                                                  Method_Mb_OrbBased_nOrbs(bt2), &
                                                  Method_Mb_OrbBased_nOrbs(bt3)))

    allocate (coeffsTmp, source=coeffs)

    do j3 = 1, size(rdm3Bt, 6)
      do j2 = 1, size(rdm3Bt, 5)
        do j1 = 1, size(rdm3Bt, 4)
          do i3 = 1, size(rdm3Bt, 3)
            do i2 = 1, size(rdm3Bt, 2)
              do i1 = 1, size(rdm3Bt, 1)

                coeffsTmp = coeffs

                call ApplyGroupedExcitations(coeffsTmp, [i1, i2, i3], [j1, j2, j3], [bt1, bt2, bt3], 3)

                rdm3Bt(i1, i2, i3, j1, j2, j3) = dot_product(coeffs, coeffsTmp)
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Compute 2-body RDM for specified body-type pair.
!>
!> Computes:
!>   ρ(i1,i2,j1,j2) = ⟨Ψ| a†_{i1,bt1} a†_{i2,bt2} a_{j2,bt2} a_{j1,bt1} |Ψ⟩
!>
!> For same body type (bt1 == bt2), this gives the intra-species 2-RDM.
!> For different body types, this gives the inter-species correlation matrix.
!>
!> **Complexity:** O(n₁² n₂² · nCoeffs)
!>
!> @param[out] rdm2Bt  4D array (n1×n2×n1×n2) allocated on output
!> @param[in]  coeffs  CI coefficient vector
!> @param[in]  bt1,bt2 Body type indices (1..nBodyTypes)
  subroutine FillRdm2Bt(rdm2Bt, coeffs, bt1, bt2)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm2Bt(:, :, :, :)
    complex(R64), intent(in), contiguous, target     :: coeffs(:)
    integer(I32), intent(in) :: bt1
    integer(I32), intent(in) :: bt2

    complex(R64), allocatable  :: coeffsTmp(:)

    integer(I32) :: i1, i2, j1, j2

    if (.not. allocated(rdm2Bt)) allocate (rdm2Bt(Method_Mb_OrbBased_nOrbs(bt1), &
                                                  Method_Mb_OrbBased_nOrbs(bt2), &
                                                  Method_Mb_OrbBased_nOrbs(bt1), &
                                                  Method_Mb_OrbBased_nOrbs(bt2)))

    allocate (coeffsTmp, source=coeffs)

    do j2 = 1, size(rdm2Bt, 4)
      do j1 = 1, size(rdm2Bt, 3)
        do i2 = 1, size(rdm2Bt, 2)
          do i1 = 1, size(rdm2Bt, 1)

            coeffsTmp = coeffs

            call ApplyGroupedExcitations(coeffsTmp, [i1, i2], [j1, j2], [bt1, bt2], 2)

            rdm2Bt(i1, i2, j1, j2) = dot_product(coeffs, coeffsTmp)
          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Compute 1-body RDM for a single body type.
!>
!> Computes:
!>   ρ(i,j) = ⟨Ψ| a†_{i,bt} a_{j,bt} |Ψ⟩
!>
!> The diagonal elements ρ(i,i) give orbital occupations; eigenvalues are
!> natural occupation numbers and eigenvectors are natural orbitals.
!>
!> **Complexity:** O(n² · nCoeffs)
!>
!> @param[out] rdm1Bt  2D array (n×n) allocated on output
!> @param[in]  coeffs  CI coefficient vector
!> @param[in]  bt1     Body type index (1..nBodyTypes)
  subroutine FillRdm1Bt(rdm1Bt, coeffs, bt1)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm1Bt(:, :)
    complex(R64), intent(in), contiguous, target     :: coeffs(:)
    integer(I32), intent(in) :: bt1

    complex(R64), allocatable  :: coeffsTmp(:)

    integer(I32) :: i1, j1

    if (.not. allocated(rdm1Bt)) allocate (rdm1Bt(Method_Mb_OrbBased_nOrbs(bt1), &
                                                  Method_Mb_OrbBased_nOrbs(bt1)))

    allocate (coeffsTmp, source=coeffs)

    do j1 = 1, size(rdm1Bt, 2)
      do i1 = 1, size(rdm1Bt, 1)

        coeffsTmp = coeffs

        call ApplyGroupedExcitations(coeffsTmp, [i1], [j1], [bt1], 1)

        rdm1Bt(i1, j1) = dot_product(coeffs, coeffsTmp)

      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
