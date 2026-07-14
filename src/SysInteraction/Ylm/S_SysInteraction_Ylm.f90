! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for Ylm-domain interactions.
!>
!> @details Provides the common spherical-harmonic routines:
!>   - `FillInteractionSrc`: Product of orbitals expanded to higher lmax
!>   - `FillInteractionPotential`: Loop over (l,m) channels, solve radial Poisson
!>   - `MultiplyWithInteractionPotential`: Angular recoupling via SpatialProduct
!>
!> The radial solver (`FillInteractionPotentialRadial`) is set by the specific
!> implementation (StdImpl, TwoScan, FullEq, BlockEq).
submodule(M_SysInteraction_Ylm) S_SysInteraction_Ylm

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Fabricate the Ylm interaction back-end.
  !>
  !> Reads `lmax` for the potential expansion (defaults to 2× orbital lmax to
  !> capture the full density product) and dispatches to Coulomb variants.
  module subroutine SysInteraction_Ylm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_SphericalHarmonics
    use M_Grid_Ylm
    use M_SysInteraction
    use M_SysInteraction_Ylm_Coulomb

    call Say_Fabricate("sysInteraction.ylm")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Ylm_lmax = Json_Get("sysInteraction.ylm.lmax", 2 * Grid_Ylm_lmax)

    ! Precompute Gaunt coefficients covering both SpatialProduct usages:
    ! source fill (l1,l2 <= lmax, l3 <= lmaxPot) and potential multiply
    ! (l1 <= lmaxPot, l2,l3 <= lmax)
    call SphericalHarmonics_EnsureGauntTable(max(Grid_Ylm_lmax, SysInteraction_Ylm_lmax), &
                                             Grid_Ylm_lmax, &
                                             max(Grid_Ylm_lmax, SysInteraction_Ylm_lmax))

    SysInteraction_FillInteractionPotential => FillInteractionPotential
    SysInteraction_FillInteractionSrc => FillInteractionSrc
    SysInteraction_MultiplyWithInteractionPotential => MultiplyWithInteractionPotential
    SysInteraction_ConjugateInteractionPotential => ConjugateInteractionPotential

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysInteraction.ylm.coulomb")) then
      call SysInteraction_Ylm_Coulomb_Fabricate

    else
      error stop "sysInteraction.ylm is missing one of: coulomb"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Compute interaction potential from source via channel-by-channel Poisson solve.
  !>
  !> For each (l,m) channel:
  !>   1. Extract radial component ρₗₘ(r) from source
  !>   2. Solve radial Poisson equation → Vₗₘ(r)
  !>   3. Store in full Ylm expansion
  !>
  !> Channels with negligible source (|ρₗₘ| < 10⁻¹⁴) are skipped for efficiency.
  !>
  !> @param[out] interactionPotential  Full Ylm-expanded potential V(r,Ω)
  !> @param[in]  src                   Source density in Ylm basis (higher lmax)
  !> @param[in]  time                  Physical time
  !> @param[in]  bt1_                  Target body type (optional)
  !> @param[in]  bt2_                  Source body type (optional)
  subroutine FillInteractionPotential(interactionPotential, src, time, bt1_, bt2_)
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), allocatable :: interactionPotential(:)
    complex(R64), intent(in), contiguous :: src(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt1_
    integer(I32), intent(in), optional :: bt2_

    integer(I32) :: l, m, potSize, nRad
    integer(I32) :: lmaxPot
    integer(I32) :: iCh, nCh
    integer(I32), allocatable :: lOfCh(:), mOfCh(:)
    complex(R64), allocatable :: srcLm(:), potLm(:)

    nRad = Grid_Ylm_nRadial
    lmaxPot = SysInteraction_Ylm_lmax
    potSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    if (.not. allocated(interactionPotential)) allocate (interactionPotential(potSize))
    interactionPotential = 0.0_R64

    ! Flatten the (l,m) channels so the independent radial Poisson solves can
    ! run in parallel; each channel writes a disjoint slice of the potential
    nCh = (lmaxPot + 1)**2
    allocate (lOfCh(nCh), mOfCh(nCh))
    iCh = 0
    do l = 0, lmaxPot
      do m = -l, l
        iCh = iCh + 1
        lOfCh(iCh) = l
        mOfCh(iCh) = m
      end do
    end do

    !$omp parallel default(shared) private(iCh, l, m, srcLm, potLm)
    allocate (srcLm(nRad), potLm(nRad))
    !$omp do schedule(dynamic)
    do iCh = 1, nCh
      l = lOfCh(iCh)
      m = mOfCh(iCh)
      call Grid_Ylm_GetLmComponent(srcLm, l, m, src)
      if (all(abs(srcLm) < 1.0e-14_R64)) cycle
      call SysInteraction_Ylm_FillInteractionPotentialRadial(potLm, srcLm, l, m, time, bt1_, bt2_)
      call Grid_Ylm_SetLmComponent(interactionPotential, l, m, potLm)
    end do
    !$omp end do
    deallocate (srcLm, potLm)
    !$omp end parallel

    deallocate (lOfCh, mOfCh)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Build source term from orbital product in Ylm basis.
  !>
  !> Computes ρ(r,Ω) = ψ*ᵢ(r,Ω) ψⱼ(r,Ω) expanded to lmax_pot to capture the
  !> full angular content of the density (product of two l functions needs 2l).
  !>
  !> Uses `Grid_Ylm_SpatialProduct` which handles the Clebsch-Gordan coupling
  !> of spherical harmonics and includes radial quadrature weights.
  !>
  !> @param[out] src        Source density ρₗₘ(r) in Ylm basis
  !> @param[in]  orbConjg   Conjugated orbital (will be conjugated → bra)
  !> @param[in]  orb        Orbital (ket)
  subroutine FillInteractionSrc(src, orbConjg, orb)
    use M_Grid_Ylm

    complex(R64), intent(out), allocatable :: src(:)
    complex(R64), intent(in), contiguous :: orbConjg(:)
    complex(R64), intent(in), contiguous :: orb(:)

    integer(I32) :: lmax, srcSize, lmaxPot

    lmax = Grid_Ylm_lmax
    lmaxPot = SysInteraction_Ylm_lmax
    srcSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    if (.not. allocated(src)) allocate (src(srcSize))

    call Grid_Ylm_SpatialProduct(src, orbConjg, orb, lmaxPot, lmax, lmax, conjgQ_=.true., withWeightsQ_=.true.)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply interaction potential to an orbital via angular recoupling.
  !>
  !> Computes dψ(r,Ω) = V(r,Ω) × ψ(r,Ω) where both are in Ylm basis.
  !> Uses `Grid_Ylm_SpatialProduct` to handle the Clebsch-Gordan coupling.
  !>
  !> @param[out] dOrb                   Resulting orbital V·ψ
  !> @param[in]  interactionPotential   Potential V in Ylm basis (lmax_pot)
  !> @param[in]  orb                    Source orbital ψ in Ylm basis (lmax)
  subroutine MultiplyWithInteractionPotential(dOrb, interactionPotential, orb)
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: interactionPotential(:)
    complex(R64), intent(in), contiguous :: orb(:)

    integer(I32) :: lmax, lmaxPot

    lmax = Grid_Ylm_lmax
    lmaxPot = SysInteraction_Ylm_lmax

    call Grid_Ylm_SpatialProduct(dOrb, interactionPotential, orb, lmax, lmaxPot, lmax)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Conjugate an interaction potential in the Ylm representation.
  !>
  !> Complex conjugation of a spatial function acts on Ylm expansion
  !> coefficients as f*_{lm} = (-1)^m conjg(f_{l,-m}), because
  !> conjg(Y_{lm}) = (-1)^m Y_{l,-m}.
  subroutine ConjugateInteractionPotential(potOut, potIn)
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous :: potOut(:)
    complex(R64), intent(in), contiguous  :: potIn(:)

    integer(I32) :: l, m, phase, nRad, lmaxPot
    complex(R64), allocatable :: potLm(:)

    nRad = Grid_Ylm_nRadial
    lmaxPot = SysInteraction_Ylm_lmax

    potOut(:) = (0.0_R64, 0.0_R64)

    allocate (potLm(nRad))

    do l = 0, lmaxPot
      do m = -l, l
        call Grid_Ylm_GetLmComponent(potLm, l, -m, potIn)
        if (mod(abs(m), 2) .eq. 0) then
          phase = 1
        else
          phase = -1
        end if
        potLm(:) = phase * conjg(potLm(:))
        call Grid_Ylm_SetLmComponent(potOut, l, m, potLm)
      end do
    end do

    deallocate (potLm)

  end subroutine

end submodule
