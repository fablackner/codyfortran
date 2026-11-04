! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Mean field of a factored one-body density matrix on top of the
!> interaction contract: direct potential, adaptively compressed exchange
!> (ACE) and their application.
!>
!> @details For gamma(r, r') = sum_l g_l(r) conjg(phi_l(r')) (phi = orbs,
!> g = weightedOrbs) the mean field is J - K with
!>   J = V[rho], rho = sum_l conjg(phi_l) g_l: SysInteraction_FillDirectPotential
!>       (one potential solve of the total source),
!>   K v = sum_l V[conjg(phi_l) v] g_l, compressed to
!>       K_ACE = sum_k |zeta_k><zeta_k|: SysInteraction_FillExchangeAce,
!> and SysInteraction_AddMeanField applies J - K_ACE to a vector. They use
!> only the SysInteraction procedure pointers and Grid_InnerProduct, so they
!> work on every grid. Users: groundSolver.scf.stdImpl (closed shell: J of
!> 2 sum |phi_l|^2 with factor_ = 2; aceQ: K_ACE of the occupied orbitals of
!> one spin rebuilt in every SCF iteration, applied with J in every Fock
!> action) and the frozen mean-field reference of the Cayley exact part of
!> Method (method.mb.orbBased.exactFieldFreeMeanFieldQ: J and K of the 1-RDM
!> of the state, built once, applied per vector).
submodule(M_SysInteraction) S_SysInteraction_MeanField

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_FillDirectPotential(pot, orbs, weightedOrbs, time, factor_, bt1_, bt2_)
    !---------------------------------------------------------------------------
    ! pot = V[src] of the total source src = factor sum_l src_l with
    ! src_l = SysInteraction_FillInteractionSrc(orbs_l, weightedOrbs_l),
    ! accumulated in the order l = 1, 2, ... (one potential solve: the
    ! potential is linear in the source). factor 1 (default) and 2 are exact
    ! scalings.
    !---------------------------------------------------------------------------
    complex(R64), intent(out), allocatable :: pot(:)
    complex(R64), intent(in), contiguous   :: orbs(:, :)
    complex(R64), intent(in), contiguous   :: weightedOrbs(:, :)
    real(R64), intent(in)                  :: time
    real(R64), intent(in), optional        :: factor_
    integer(I32), intent(in), optional     :: bt1_
    integer(I32), intent(in), optional     :: bt2_

    complex(R64), allocatable :: src(:), srcTotal(:)
    real(R64) :: factor
    integer(I32) :: l

    if (any(shape(weightedOrbs) .ne. shape(orbs))) &
      error stop "SysInteraction_FillDirectPotential: orbs and weightedOrbs must have the same shape"
    if (size(orbs, 2) < 1) error stop "SysInteraction_FillDirectPotential: at least one orbital is needed"

    factor = 1.0_R64
    if (present(factor_)) factor = factor_

    do l = 1, size(orbs, 2)
      call SysInteraction_FillInteractionSrc(src, orbs(:, l), weightedOrbs(:, l))
      if (l .eq. 1) then
        srcTotal = factor * src
      else
        srcTotal = srcTotal + factor * src
      end if
    end do
    call SysInteraction_FillInteractionPotential(pot, srcTotal, time, bt1_=bt1_, bt2_=bt2_)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_FillExchangeAce(zeta, positiveDefiniteQ, orbs, weightedOrbs, time, bt1_, bt2_)
    !---------------------------------------------------------------------------
    ! For the n orbitals phi_i = orbs(:, i) and the weighted orbitals
    ! g_l = weightedOrbs(:, l):
    !   xi_j = K phi_j = sum_k V_kj g_k,  V_kj = V[conjg(phi_k) phi_j],
    ! with V_kj solved for k <= j only when SysInteraction_conjSymmetricQ
    ! (V_jk = conjg(V_kj) via SysInteraction_ConjugateInteractionPotential),
    ! else for all pairs;
    !   M_ij = <phi_i|xi_j> (Hermitized), M = L L^H (Cholesky), zeta = xi L^-H,
    ! so that sum_k zeta_k <zeta_k|phi> = xi M^-1 <xi|phi> reproduces K phi for
    ! phi in span{phi_j}. If M is numerically not positive definite (e.g. a
    ! slightly non-symmetric discretized interaction kernel), the eigen-
    ! decomposition M = U Lambda U^H is used instead: zeta = xi U Lambda^(-1/2)
    ! on the eigenvalues lambda > 1e-12 lambda_max (the other columns of zeta
    ! are zero), i.e. the pseudo-inverse of M, and positiveDefiniteQ is false.
    !---------------------------------------------------------------------------
    use M_Utils_LapackLib, only: ZPOTRF, LapackLib_DiagonalizeGeneric
    use M_Utils_BlasLib, only: ZGEMM
    use M_Grid

    complex(R64), intent(out), contiguous :: zeta(:, :)
    logical, intent(out)                  :: positiveDefiniteQ
    complex(R64), intent(in), contiguous  :: orbs(:, :)
    complex(R64), intent(in), contiguous  :: weightedOrbs(:, :)
    real(R64), intent(in)                 :: time
    integer(I32), intent(in), optional    :: bt1_
    integer(I32), intent(in), optional    :: bt2_

    complex(R64), allocatable :: xi(:, :), mat(:, :), matCopy(:, :), lInv(:, :), u(:, :)
    complex(R64), allocatable :: src(:), pot(:), potConjg(:), tmp(:)
    real(R64), allocatable :: lambda(:)
    integer(I32) :: nG, n, i, j, k, kMax, info

    nG = int(size(orbs, 1), I32)
    n = int(size(orbs, 2), I32)

    if (any(shape(weightedOrbs) .ne. shape(orbs)) .or. any(shape(zeta) .ne. shape(orbs))) &
      error stop "SysInteraction_FillExchangeAce: orbs, weightedOrbs and zeta must have the same shape"

    positiveDefiniteQ = .true.
    if (n .eq. 0) return

    ! ξ_j = K̂φ_j
    allocate (xi(nG, n), tmp(nG))
    xi(:, :) = (0.0_R64, 0.0_R64)
    do j = 1, n
      if (SysInteraction_conjSymmetricQ) then
        kMax = j
      else
        kMax = n
      end if
      do k = 1, kMax
        ! V_kj = V[conjg(φ_k) φ_j] contributes V_kj·g_k to ξ_j
        call SysInteraction_FillInteractionSrc(src, orbs(:, k), orbs(:, j))
        call SysInteraction_FillInteractionPotential(pot, src, time, bt1_=bt1_, bt2_=bt2_)
        call SysInteraction_MultiplyWithInteractionPotential(tmp, pot, weightedOrbs(:, k))
        xi(:, j) = xi(:, j) + tmp
        if (SysInteraction_conjSymmetricQ .and. k .ne. j) then
          ! V_jk = conj(V_kj) contributes V_jk·g_j to ξ_k
          if (allocated(potConjg)) then
            if (size(potConjg) .ne. size(pot)) deallocate (potConjg)
          end if
          if (.not. allocated(potConjg)) allocate (potConjg, mold=pot)
          call SysInteraction_ConjugateInteractionPotential(potConjg, pot)
          call SysInteraction_MultiplyWithInteractionPotential(tmp, potConjg, weightedOrbs(:, j))
          xi(:, k) = xi(:, k) + tmp
        end if
      end do
    end do

    ! M_ij = ⟨φ_i|ξ_j⟩, Hermitian part
    allocate (mat(n, n))
    do j = 1, n
      do i = 1, n
        mat(i, j) = Grid_InnerProduct(orbs(:, i), xi(:, j))
      end do
    end do
    mat = 0.5_R64 * (mat + conjg(transpose(mat)))
    matCopy = mat

    ! M = L Lᴴ and L⁻¹ (forward substitution per column)
    info = 0
    call ZPOTRF('L', n, mat, n, info)
    allocate (lInv(n, n))
    lInv(:, :) = (0.0_R64, 0.0_R64)
    positiveDefiniteQ = info .eq. 0

    if (positiveDefiniteQ) then
      do j = 1, n
        lInv(j, j) = 1.0_R64 / mat(j, j)
        do i = j + 1, n
          lInv(i, j) = -sum(mat(i, j:i - 1) * lInv(j:i - 1, j)) / mat(i, i)
        end do
      end do

      ! ζ = ξ L⁻ᴴ
      call ZGEMM('N', 'C', nG, n, n, (1.0_R64, 0.0_R64), xi, nG, lInv, n, (0.0_R64, 0.0_R64), zeta, nG)
    else
      ! Not numerically positive definite: ζ = ξ U Λ^(-1/2) on the positive part
      call LapackLib_DiagonalizeGeneric(lambda, u, matCopy, .true.)
      do j = 1, n
        if (lambda(j) > 1e-12_R64 * maxval(lambda)) lInv(:, j) = u(:, j) / sqrt(lambda(j))
      end do
      call ZGEMM('N', 'N', nG, n, n, (1.0_R64, 0.0_R64), xi, nG, lInv, n, (0.0_R64, 0.0_R64), zeta, nG)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_AddMeanField(w, scratch, z, directPotential, zeta)
    !---------------------------------------------------------------------------
    ! w = w + J z - sum_k zeta_k <zeta_k|z>: the overlaps first (in parallel
    ! over k), J z in scratch, then one fused pass, per element
    ! ((w_i + (J z)_i) - o_1 zeta_i1) - o_2 zeta_i2 - ... Reentrant; its
    ! passes are parallel when called outside a parallel region.
    !---------------------------------------------------------------------------
    !$ use omp_lib, only: omp_in_parallel
    use M_Grid

    complex(R64), intent(inout), contiguous :: w(:)
    complex(R64), intent(inout), contiguous :: scratch(:)
    complex(R64), intent(in), contiguous    :: z(:)
    complex(R64), intent(in), contiguous    :: directPotential(:)
    complex(R64), intent(in), contiguous    :: zeta(:, :)

    complex(R64), allocatable :: overlaps(:)
    complex(R64) :: acc
    integer(I32) :: i, k, nZeta

    if (size(z) .ne. size(w) .or. size(scratch) .ne. size(w) .or. size(zeta, 1) .ne. size(w)) &
      error stop "SysInteraction_AddMeanField: w, scratch, z and the columns of zeta must have the same size"

    nZeta = size(zeta, 2)
    allocate (overlaps(nZeta))
    ! The inner products in parallel over k (Grid_InnerProduct is reentrant
    ! and, e.g. on grid.fem, a serial pass)
    !$omp parallel do default(shared) private(k) schedule(static) if (.not. omp_in_parallel())
    do k = 1, nZeta
      overlaps(k) = Grid_InnerProduct(zeta(:, k), z)
    end do
    !$omp end parallel do

    call SysInteraction_MultiplyWithInteractionPotential(scratch, directPotential, z)

    !$omp parallel do default(shared) private(i, k, acc) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(w)
      acc = w(i) + scratch(i)
      do k = 1, nZeta
        acc = acc - overlaps(k) * zeta(i, k)
      end do
      w(i) = acc
    end do
    !$omp end parallel do

  end subroutine

end submodule
