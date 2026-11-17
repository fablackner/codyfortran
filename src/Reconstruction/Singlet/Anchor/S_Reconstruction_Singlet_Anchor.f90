! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The anchored closure F(D) = ccQ(R(D) + T_S delta3), delta3 = D3_ref -
!> R(D2_ref), and the spectral part; see A_Reconstruction_Singlet_Anchor.md.
submodule(M_Reconstruction_Singlet_Anchor) S_Reconstruction_Singlet_Anchor
  use M_Reconstruction, only: I_Reconstruction_AddCollisionOp
  use M_Reconstruction_Singlet, only: I_Reconstruction_Singlet_FillRdm3aab

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Reference files (in anchorDir): the 2-RDM of Coeffs_SaveTwoRdm, the
  !> 3-RDM and transition RDMs of Reconstruction_Singlet_Anchor_SaveReference.
  character(len=*), parameter :: rdm2File = "rdm2Bt01_02.in", rdm3File = "rdm3Bt01_01_02.in"
  character(len=*), parameter :: trdm2File = "trdm2Bt01_02.in", trdm3File = "trdm3Bt01_01_02.in"

  !> The closure of the scheme and the collision operator of the family that
  !> this extension wraps.
  procedure(I_Reconstruction_Singlet_FillRdm3aab), pointer :: SchemeFillRdm3aab => null()
  procedure(I_Reconstruction_AddCollisionOp), pointer :: FamilyAddCollisionOp => null()

  !> Directory of the reference files, "" or with a trailing "/".
  character(len=:), allocatable :: anchorDir
  !> Number of transitions of the spectral part (0: reference only).
  integer(I32) :: nTransitions = 0
  !> Orbitals of the reference (Orbs_orbs at Setup) and the overlap
  !> S(i, j) = <phi_i|phi_ref,j> of the current ones (Reconstruction_SetOrbitals;
  !> unallocated for static orbitals: identity).
  complex(R64), allocatable :: refOrbs(:, :), overlap(:, :)
  !> The 2-RDM of the reference and delta3 in the orbitals of the reference.
  complex(R64), allocatable :: anchorRdm2ab(:, :, :, :), anchorDelta3(:, :, :, :, :, :)
  !> Spectral part (static orbitals): conjugated duals conjg(u_a) (n^4 x 2K),
  !> 3-RDM corrections e_a (n^6 x 2K), their collision operators (n^4 x 2K)
  !> for the h2 in spectralH2; the correction of a 2-RDM D is
  !> sum_a Re<u_a, D - D_ref> e_a.
  complex(R64), allocatable :: spectralDuals(:, :), spectralCorrections(:, :)
  complex(R64), allocatable :: spectralCollisions(:, :), spectralH2(:, :, :, :)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Singlet_Anchor_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Reconstruction
    use M_Reconstruction_Singlet

    call Say_Fabricate("reconstruction.singlet.anchor")

    anchorDir = Json_Get("reconstruction.singlet.anchor.dir", "")
    if (len_trim(anchorDir) > 0) anchorDir = trim(anchorDir)//"/"
    nTransitions = Json_Get("reconstruction.singlet.anchor.nTransitions", 0)
    if (nTransitions < 0) error stop "reconstruction.singlet.anchor.nTransitions must be >= 0"
    call ClearState

    SchemeFillRdm3aab => Reconstruction_Singlet_FillRdm3aab
    Reconstruction_Singlet_FillRdm3aab => FillRdm3aab
    if (nTransitions > 0) then
      FamilyAddCollisionOp => Reconstruction_AddCollisionOp
      Reconstruction_AddCollisionOp => AddCollisionOp
    end if
    Reconstruction_Setup => Setup
    Reconstruction_SetOrbitals => SetOrbitals

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Singlet_Anchor_SaveReference(kets)
    use M_Utils_DataStorage
    use M_Coeffs

    complex(R64), intent(in), contiguous :: kets(:, :)

    complex(R64), allocatable :: rdm2Bt(:, :, :, :, :), rdm3Bt(:, :, :, :, :, :, :)

    call Coeffs_FillTransitionRdm3Bt(rdm3Bt, kets(:, 1), kets, 1, 1, 2)
    call SaveData(rdm3Bt(:, :, :, :, :, :, 1), rdm3File, storage_size(rdm3Bt))
    if (size(kets, 2) .eq. 1) return
    call SaveData(rdm3Bt(:, :, :, :, :, :, 2:), trdm3File, storage_size(rdm3Bt))
    call Coeffs_FillTransitionRdm2Bt(rdm2Bt, kets(:, 1), kets(:, 2:), 1, 2)
    call SaveData(rdm2Bt, trdm2File, storage_size(rdm2Bt))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ClearState

    if (allocated(refOrbs)) deallocate (refOrbs)
    if (allocated(overlap)) deallocate (overlap)
    if (allocated(anchorRdm2ab)) deallocate (anchorRdm2ab)
    if (allocated(anchorDelta3)) deallocate (anchorDelta3)
    if (allocated(spectralDuals)) deallocate (spectralDuals)
    if (allocated(spectralCorrections)) deallocate (spectralCorrections)
    if (allocated(spectralCollisions)) deallocate (spectralCollisions)
    if (allocated(spectralH2)) deallocate (spectralH2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Reads the reference RDMs, stores delta3 = D3_ref - R(D2_ref) and the
  !> reference orbitals; with nTransitions > 0 also the spectral part.
  subroutine Setup
    use M_Utils_Say
    use M_Utils_DataStorage
    use M_Utils_RdmTrace
    use M_Orbs
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), allocatable :: rdm1a(:, :), rdm3aab(:, :, :, :, :, :)
    real(R64) :: refNorm
    integer(I32) :: n

    call Say_Setup("reconstruction.singlet.anchor")
    call ClearState

    n = Method_Mb_OrbBased_nOrbs(1)
    refOrbs = Orbs_orbs
    allocate (anchorRdm2ab(n, n, n, n), anchorDelta3(n, n, n, n, n, n))
    call LoadData(anchorRdm2ab, anchorDir//rdm2File, storage_size(anchorRdm2ab))
    call LoadData(anchorDelta3, anchorDir//rdm3File, storage_size(anchorDelta3))

    call RdmTrace_FillRdm1aFromRdm2ab(rdm1a, anchorRdm2ab, Method_Mb_nBodies(2))
    call SchemeFillRdm3aab(rdm3aab, rdm1a, anchorRdm2ab)
    refNorm = norm2(abs(anchorDelta3))
    anchorDelta3(:, :, :, :, :, :) = anchorDelta3(:, :, :, :, :, :) - rdm3aab(:, :, :, :, :, :)
    write (*, '(A,ES12.4,A,ES12.4)') "  anchor: |D3_ref - R(D2_ref)| = ", norm2(abs(anchorDelta3)), &
      ",  |D3_ref| = ", refNorm

    if (nTransitions > 0) call SetupSpectral(n)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Spectral part from the transition RDMs T2_k = <ref|..|chi_k>, T3_k
  !> (k = 1..nTransitions): the Hermitian perturbations v_a = t + t^+ of the
  !> 2-RDM, t = T2_k and i T2_k, and e_a = (exact 3-RDM change) - (change of
  !> F, central difference); the minimal-norm correction exact on span{v_a}
  !> has the duals u_a = sum_b G+_ab v_b, G_ab = Re<v_a, v_b>.
  subroutine SetupSpectral(n)
    use M_Utils_Constants, only: CONE, IU
    use M_Utils_DataStorage
    use M_Utils_LapackLib
    use M_Utils_RdmEnforceTraceM3aab
    use M_Reconstruction_Singlet

    integer(I32), intent(in) :: n

    complex(R64), allocatable :: trdm2(:, :, :, :, :), trdm3(:, :, :, :, :, :), v(:, :)
    complex(R64), allocatable :: rdm3Plus(:, :, :, :, :, :), rdm3Minus(:, :, :, :, :, :)
    real(R64), allocatable :: evals(:), evecs(:, :)
    real(R64), parameter :: eps = 1e-4_R64
    complex(R64) :: phase
    integer(I32) :: k, b, a, nV

    nV = 2 * nTransitions
    allocate (trdm2(n, n, n, n, nTransitions), trdm3(n, n, n, n, n, n), v(n**4, nV), spectralCorrections(n**6, nV))
    call LoadData(trdm2, anchorDir//trdm2File, storage_size(trdm2))
    do k = 1, nTransitions
      call LoadData(trdm3, anchorDir//trdm3File, storage_size(trdm3), startPos_=(k - 1) * n**6 + 1)
      do b = 1, 2
        a = 2 * (k - 1) + b
        phase = merge(CONE, IU, b .eq. 1)
        v(:, a) = reshape(phase * trdm2(:, :, :, :, k) + &
                          conjg(phase * reshape(trdm2(:, :, :, :, k), [n, n, n, n], order=[3, 4, 1, 2])), [n**4])
        call Reconstruction_Singlet_ReconstructRdm3aab(rdm3Plus, anchorRdm2ab + eps * reshape(v(:, a), [n, n, n, n]))
        call Reconstruction_Singlet_ReconstructRdm3aab(rdm3Minus, anchorRdm2ab - eps * reshape(v(:, a), [n, n, n, n]))
        rdm3Plus = phase * trdm3 + conjg(phase * reshape(trdm3, [n, n, n, n, n, n], order=[4, 5, 6, 1, 2, 3])) - &
                   (rdm3Plus - rdm3Minus) / (2.0_R64 * eps)
        ! with ccQ, e_a has zero contractions up to the difference error (F and
        ! the exact change both carry those of v_a): project exactly, so that it
        ! may enter after ccQ; without ccQ it enters as it is
        if (Reconstruction_Singlet_ccQ) call RdmEnforceTraceM3aab_Apply(rdm3Plus)
        spectralCorrections(:, a) = reshape(rdm3Plus, [n**6])
      end do
    end do

    call LapackLib_DiagonalizeGeneric(evals, evecs, real(matmul(conjg(transpose(v)), v), R64), .true.)
    write (*, '(A,I0,A,ES12.4,A,ES12.4)') "  spectral anchor: ", nTransitions, " transitions,  max |e_a| = ", &
      maxval(norm2(abs(spectralCorrections), dim=1)), ",  min/max eig of the Gram matrix ", minval(evals) / maxval(evals)
    where (evals > 1e-10_R64 * maxval(evals))
      evals = 1.0_R64 / evals
    elsewhere
      evals = 0.0_R64
    end where
    spectralDuals = conjg(matmul(v, cmplx(matmul(evecs, spread(evals, 2, nV) * transpose(evecs)), 0.0_R64, R64)))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> The closure of the scheme plus delta3 in the current orbitals (S on the
  !> ket, conjg(S) on the bra indices; identity for static orbitals).
  subroutine FillRdm3aab(rdm3aab, rdm1a, rdm2ab)
    use M_Utils_RdmBasisChange

    complex(R64), intent(out), allocatable :: rdm3aab(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm1a(:, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    complex(R64), allocatable :: delta3(:, :, :, :, :, :)

    if (.not. allocated(anchorDelta3)) &
      error stop "reconstruction.singlet.anchor: call Reconstruction_Setup before Reconstruction_AddCollisionOp"

    call SchemeFillRdm3aab(rdm3aab, rdm1a, rdm2ab)
    if (allocated(overlap)) then
      delta3 = anchorDelta3
      call RdmBasisChange_Transform(delta3, 6, overlap)
      rdm3aab(:, :, :, :, :, :) = rdm3aab(:, :, :, :, :, :) + delta3(:, :, :, :, :, :)
    else
      rdm3aab(:, :, :, :, :, :) = rdm3aab(:, :, :, :, :, :) + anchorDelta3(:, :, :, :, :, :)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> The overlap of the current orbitals with the reference ones (moving
  !> orbitals); the spectral part needs static orbitals.
  subroutine SetOrbitals(orbs)
    use M_Method_Mb_OrbBased

    complex(R64), intent(in), contiguous :: orbs(:, :)

    complex(R64), allocatable :: s(:, :)

    if (nTransitions > 0) error stop "reconstruction.singlet.anchor.nTransitions > 0 needs static orbitals (td2rdmStaticOrbs)"
    if (.not. allocated(refOrbs)) error stop "reconstruction.singlet.anchor: call Reconstruction_Setup first"

    call Method_Mb_OrbBased_FillOneBodyMatrix(s, orbs, refOrbs)
    overlap = s(1:Method_Mb_OrbBased_nOrbs(1), 1:Method_Mb_OrbBased_nOrbs(1))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> The collision operator of the family plus sum_a c_a Tr3[W, e_a],
  !> c_a = Re<u_a, D - D_ref>: ccQ is affine and leaves the e_a unchanged, the
  !> collision operator is real-linear, so the spectral part enters through
  !> the collision operators of the e_a, computed whenever h2 changes.
  subroutine AddCollisionOp(dTwoRdmState, twoRdmState, h2)
    use M_Utils_RdmCollisionOp
    use M_TwoRdm_Fermionic_Singlet

    complex(R64), intent(inout), contiguous :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous    :: twoRdmState(:)
    complex(R64), intent(in), contiguous    :: h2(:, :, :, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :), dRdm2ab(:, :, :, :), collOp2ab(:, :, :, :)
    real(R64), allocatable :: coefficients(:)
    logical :: cachedQ
    integer(I32) :: n, a

    call FamilyAddCollisionOp(dTwoRdmState, twoRdmState, h2)

    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(dRdm2ab, dTwoRdmState)
    n = size(rdm2ab, 1)

    cachedQ = .false.
    if (allocated(spectralH2)) then
      if (all(shape(spectralH2) .eq. shape(h2))) cachedQ = maxval(abs(spectralH2 - h2)) <= 0.0_R64
    end if
    if (.not. cachedQ) then
      spectralH2 = h2
      if (allocated(spectralCollisions)) deallocate (spectralCollisions)
      allocate (spectralCollisions(n**4, size(spectralCorrections, 2)))
      do a = 1, size(spectralCorrections, 2)
        call RdmCollisionOp_FillCollisionOp2ab(collOp2ab, h2, reshape(spectralCorrections(:, a), [n, n, n, n, n, n]))
        spectralCollisions(:, a) = reshape(collOp2ab, [n**4])
      end do
    end if

    coefficients = real(matmul(reshape(rdm2ab - anchorRdm2ab, [n**4]), spectralDuals), R64)
    dRdm2ab(:, :, :, :) = dRdm2ab(:, :, :, :) + &
                          reshape(matmul(spectralCollisions, cmplx(coefficients, 0.0_R64, R64)), [n, n, n, n])

    call TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab(dTwoRdmState, dRdm2ab)

  end subroutine

end submodule
