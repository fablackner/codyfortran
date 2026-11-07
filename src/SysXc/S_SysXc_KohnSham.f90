! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Grid-agnostic Kohn-Sham drivers of SysXc (on the contract only).
!>
!> @details
!> - SysXc_FillXcPotential: the discrete functional E_xc = Σ_p w_p e(ρ_p, σ_p)
!>   with σ from SysXc_ApplyGradient (GGA) and its exact derivative
!>     v_s = ∂e/∂ρ_s + W⁻¹GᵀW (2 ∂e/∂σ_ss Gρ_s + ∂e/∂σ_ud Gρ_s')
!>   with SysXc_AddGradientAdjoint;
!> - SysXc_FillHartreePotential: J from one SysInteraction solve of the
!>   source of the total density, E_H = ½ Σ w ρ J;
!> - SysXc_FillKohnShamPotential: v_KS,s = J + v_xc,s;
!> - SysXc_Integrate: Σ w f.
!> The point loops are OpenMP-parallel unless called in an active parallel
!> region.
submodule(M_SysXc) S_SysXc_KohnSham

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysXc_FillXcPotential(xcPot, xcEnergy, density)
    !$ use omp_lib, only: omp_in_parallel

    real(R64), intent(out), contiguous :: xcPot(:, :)
    real(R64), intent(out) :: xcEnergy
    real(R64), intent(in), contiguous :: density(:, :)

    real(R64), allocatable :: rho(:, :), sigma(:, :), vSigma(:, :), energyDensity(:)
    real(R64), allocatable :: gradRho(:, :, :), flux(:, :)
    integer(I32) :: nP, nSpin, nSigma, nDim, s, sOther, c, p, iSs, iSo

    nP = SysXc_nPoints
    nSpin = size(density, 2)
    nDim = SysXc_dimension

    if (size(density, 1) .ne. nP .or. size(xcPot, 1) .ne. nP .or. size(xcPot, 2) .ne. nSpin) &
      error stop "SysXc_FillXcPotential: arrays must be (SysXc_nPoints, nSpin)"
    if (nSpin .ne. 1 .and. nSpin .ne. 2) error stop "SysXc_FillXcPotential: nSpin must be 1 or 2"

    ! Round-off can leave tiny negative densities (weighted orbitals)
    allocate (rho(nP, nSpin))
    !$omp parallel do default(shared) private(p, s) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do p = 1, nP
      do s = 1, nSpin
        rho(p, s) = max(density(p, s), 0.0_R64)
      end do
    end do
    !$omp end parallel do

    !------------------------------------
    ! gradients and their invariants
    !------------------------------------

    nSigma = 0
    if (SysXc_gradientQ) then
      nSigma = 1
      if (nSpin .eq. 2) nSigma = 3
      allocate (gradRho(nP, nDim, nSpin))
      do s = 1, nSpin
        call SysXc_ApplyGradient(gradRho(:, :, s), rho(:, s))
      end do

      allocate (sigma(nP, nSigma))
      !$omp parallel do default(shared) private(p) schedule(static) &
      !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
      do p = 1, nP
        sigma(p, 1) = sum(gradRho(p, :, 1)**2)
        if (nSpin .eq. 2) then
          sigma(p, 2) = sum(gradRho(p, :, 1) * gradRho(p, :, 2))
          sigma(p, 3) = sum(gradRho(p, :, 2)**2)
        end if
      end do
      !$omp end parallel do
    else
      allocate (sigma(nP, 0))
    end if

    allocate (vSigma(nP, nSigma), energyDensity(nP))

    !------------------------------------
    ! functional on the points
    !------------------------------------

    call SysXc_EvaluateFunctional(energyDensity, xcPot, vSigma, rho, sigma)
    xcEnergy = SysXc_Integrate(energyDensity)

    !------------------------------------
    ! gradient terms: v_s += W^-1 G^T W (2 de/dsigma_ss G rho_s + de/dsigma_ud G rho_s')
    !------------------------------------

    if (.not. SysXc_gradientQ) return

    allocate (flux(nP, nDim))
    do s = 1, nSpin
      sOther = 3 - s
      ! sigma index of the same-spin invariant: 1 (unpolarized, up-up) or 3 (down-down)
      iSs = 1
      if (s .eq. 2) iSs = 3
      iSo = 2
      !$omp parallel do default(shared) private(p, c) schedule(static) &
      !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
      do p = 1, nP
        do c = 1, nDim
          flux(p, c) = 2.0_R64 * vSigma(p, iSs) * gradRho(p, c, s)
          if (nSpin .eq. 2) flux(p, c) = flux(p, c) + vSigma(p, iSo) * gradRho(p, c, sOther)
        end do
      end do
      !$omp end parallel do
      call SysXc_AddGradientAdjoint(xcPot(:, s), flux)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysXc_FillHartreePotential(hartreePot, hartreeEnergy, density, time)
    use M_SysInteraction

    real(R64), intent(out), contiguous :: hartreePot(:)
    real(R64), intent(out) :: hartreeEnergy
    real(R64), intent(in), contiguous :: density(:)
    real(R64), intent(in) :: time

    complex(R64), allocatable :: src(:), interactionPotential(:)

    hartreePot(:) = 0.0_R64
    hartreeEnergy = 0.0_R64
    if (.not. SysXc_hartreeQ) return

    ! The Hartree potential of the electrons needs one interaction for all
    ! spins (it is linear in the source: one solve of the total density)
    if (.not. SysInteraction_bodyTypeIndependentQ) &
      error stop "sysXc: the Hartree potential requires a body-type-independent sysInteraction"

    call SysXc_FillInteractionSrc(src, density)
    call SysInteraction_FillInteractionPotential(interactionPotential, src, time)
    call SysXc_AddInteractionPotential(hartreePot, interactionPotential)

    hartreeEnergy = 0.5_R64 * SysXc_Integrate(density * hartreePot)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysXc_FillKohnShamPotential(kohnShamPot, hartreeEnergy, xcEnergy, density, time)
    !$ use omp_lib, only: omp_in_parallel

    real(R64), intent(out), contiguous :: kohnShamPot(:, :)
    real(R64), intent(out) :: hartreeEnergy
    real(R64), intent(out) :: xcEnergy
    real(R64), intent(in), contiguous :: density(:, :)
    real(R64), intent(in) :: time

    real(R64), allocatable :: hartreePot(:)
    integer(I32) :: s, p

    allocate (hartreePot(SysXc_nPoints))
    call SysXc_FillHartreePotential(hartreePot, hartreeEnergy, sum(density, dim=2), time)
    call SysXc_FillXcPotential(kohnShamPot, xcEnergy, density)

    !$omp parallel do default(shared) private(p, s) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do p = 1, SysXc_nPoints
      do s = 1, size(kohnShamPot, 2)
        kohnShamPot(p, s) = kohnShamPot(p, s) + hartreePot(p)
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function SysXc_Integrate(f) result(res)
    real(R64), intent(in), contiguous :: f(:)
    real(R64) :: res

    if (size(f) .ne. SysXc_nPoints) error stop "SysXc_Integrate: field must have SysXc_nPoints entries"
    if (.not. allocated(SysXc_weights)) error stop "SysXc_Integrate: SysXc_Setup was not called"
    res = sum(SysXc_weights(:) * f(:))

  end function

end submodule
