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
!>   - `ApplyGroupedExcitations`: Helper to batch excitations by body type
!>
!> The excitation-based RDM blocks (Coeffs_FillTransitionRdm{2,3}Bt) come
!> from one worker, FillTransitionBlock: one excitation of the bra per element
!> and one matrix-vector product with the kets. The n^(2R) excitations of a CI
!> vector make it brute force: for small systems and validation.
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
!> (normalization, projection, I/O) are assigned unconditionally
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
    Coeffs_ProjectOnSector => NoOpProcedures_ProjectCoeffs
    Coeffs_SaveCoeffs => SaveCoeffs
    Coeffs_SaveTwoRdm => SaveTwoRdm

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
!> (written by Method_Mb_OrbBased_SaveRdm2Blocks).
  subroutine SaveTwoRdm(coeffs)
    use M_Method_Mb_OrbBased

    complex(R64), intent(in), contiguous  :: coeffs(:)

    complex(R64), allocatable :: rdm2(:, :, :, :)

    call Coeffs_ApplyH2FillRdm2(coeffs, rdm2_=rdm2)
    call Method_Mb_OrbBased_SaveRdm2Blocks(rdm2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Coeffs_FillTransitionRdm2Bt(rdm2Bt, bra, kets, bt1, bt2)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm2Bt(:, :, :, :, :)
    complex(R64), intent(in), contiguous   :: bra(:)
    complex(R64), intent(in), contiguous   :: kets(:, :)
    integer(I32), intent(in) :: bt1
    integer(I32), intent(in) :: bt2

    associate (n1 => Method_Mb_OrbBased_nOrbs(bt1), n2 => Method_Mb_OrbBased_nOrbs(bt2))
      allocate (rdm2Bt(n1, n2, n1, n2, size(kets, 2)))
    end associate
    call FillTransitionBlock(rdm2Bt, bra, kets, [bt1, bt2])

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Coeffs_FillTransitionRdm3Bt(rdm3Bt, bra, kets, bt1, bt2, bt3)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm3Bt(:, :, :, :, :, :, :)
    complex(R64), intent(in), contiguous   :: bra(:)
    complex(R64), intent(in), contiguous   :: kets(:, :)
    integer(I32), intent(in) :: bt1
    integer(I32), intent(in) :: bt2
    integer(I32), intent(in) :: bt3

    associate (n1 => Method_Mb_OrbBased_nOrbs(bt1), n2 => Method_Mb_OrbBased_nOrbs(bt2), &
               n3 => Method_Mb_OrbBased_nOrbs(bt3))
      allocate (rdm3Bt(n1, n2, n3, n1, n2, n3, size(kets, 2)))
    end associate
    call FillTransitionBlock(rdm3Bt, bra, kets, [bt1, bt2, bt3])

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> rdm(l + (k - 1) nElements) = <bra| a+_j1..a+_jR a_iR..a_i1 |kets(:, k)>,
  !> l the column-major index of (i1..iR, j1..jR), i_r and j_r orbitals of body
  !> type bts(r) (ket-first, the layout of TwoRdm): one excitation of the bra
  !> per element and one matrix-vector product with the kets, OpenMP over the
  !> elements (Coeffs_ApplyExcitation must be thread-safe).
  subroutine FillTransitionBlock(rdm, bra, kets, bts)
    use M_Utils_BlasLib
    use M_Method_Mb_OrbBased

    complex(R64), intent(out)            :: rdm(*)
    complex(R64), intent(in), contiguous :: bra(:)
    complex(R64), intent(in), contiguous :: kets(:, :)
    integer(I32), intent(in)             :: bts(:)

    complex(R64), allocatable :: w(:)
    complex(R64) :: y(size(kets, 2))
    integer(I32) :: dims(2 * size(bts)), idx(2 * size(bts)), r, nElements, l, m, rest

    r = size(bts)
    dims = [Method_Mb_OrbBased_nOrbs(bts), Method_Mb_OrbBased_nOrbs(bts)]
    nElements = product(dims)

    ! w = a+_i1..a+_iR a_jR..a_j1 |bra>, <w|ket> = <bra| a+_j1..a+_jR a_iR..a_i1 |ket>
    !$omp parallel do private(w, y, idx, rest, m) schedule(dynamic)
    do l = 1, nElements
      rest = l - 1
      do m = 1, 2 * r
        idx(m) = mod(rest, dims(m)) + 1
        rest = rest / dims(m)
      end do
      w = bra
      call ApplyGroupedExcitations(w, idx(1:r), idx(r + 1:), bts, r)
      call BlasLib_CalcMatOnVec(y, kets, w, 'C')
      rdm(l:l + (size(kets, 2) - 1) * nElements:nElements) = conjg(y)
    end do
    !$omp end parallel do

  end subroutine

end submodule
