! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Dispatch submodule of the exchange-correlation functional.
!>
!> @details
!> Binds the defaults of the nodal families, whose xc points are the grid
!> nodes (grid.fem, grid.linear: orbitals are nodal values, the metric is
!> diagonal with the weights `SysXc_weights`, the interaction source of a
!> pair density is its nodal value times the weight), then dispatches to the
!> grid family:
!>
!>   - `sysXc.fem`:    3D finite elements, LDA and GGA (nodal gradient)
!>   - `sysXc.ylm`:    radial grid × spherical harmonics (own angular
!>                     quadrature; overrides every map), LDA and GGA
!>   - `sysXc.linear`: 1D grids, LDA (1D functionals of Libxc)
!>
!> Nodal defaults:
!>   - AddDensity:     density_p += Re Σ_l conj(φ_l,p) g_l,p
!>   - Multiply:       (v·orb)_p = v_p orb_p
!>   - Interaction source / potential: src_p = ρ_p w_p, v_p += Re V_p
submodule(M_SysXc) S_SysXc

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysXc_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysXc_Fem
    use M_SysXc_Ylm
    use M_SysXc_Linear

    call Say_Fabricate("sysXc")

    !------------------------------------
    ! defaults (nodal families)
    !------------------------------------

    SysXc_hartreeQ = Json_GetExistence("sysInteraction")
    SysXc_gradientQ = .false.
    SysXc_exactExchangeFraction = 0.0_R64

    SysXc_Setup => NoOpProcedures_Setup
    SysXc_AddDensity => AddDensityNodal
    SysXc_MultiplyWithPotential => MultiplyWithPotentialNodal
    SysXc_FillInteractionSrc => FillInteractionSrcNodal
    SysXc_AddInteractionPotential => AddInteractionPotentialNodal
    nullify (SysXc_ApplyGradient, SysXc_AddGradientAdjoint)

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysXc.fem")) then
      call SysXc_Fem_Fabricate

    else if (Json_GetExistence("sysXc.ylm")) then
      call SysXc_Ylm_Fabricate

    else if (Json_GetExistence("sysXc.linear")) then
      call SysXc_Linear_Fabricate

    else
      error stop "sysXc is missing one of: fem, ylm, linear"
    end if

    if (SysXc_gradientQ .and. .not. (associated(SysXc_ApplyGradient) .and. associated(SysXc_AddGradientAdjoint))) &
      error stop "sysXc: the functional needs density gradients (GGA), which this grid family does not provide"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> density_p += Re Σ_l conj(orbs(p, l)) weightedOrbs(p, l) on the nodes.
  subroutine AddDensityNodal(density, orbs, weightedOrbs)
    !$ use omp_lib, only: omp_in_parallel

    real(R64), intent(inout), contiguous, target :: density(:)
    complex(R64), intent(in), contiguous, target :: orbs(:, :)
    complex(R64), intent(in), contiguous, target :: weightedOrbs(:, :)

    integer(I32) :: p, l
    real(R64) :: acc

    if (size(density) .ne. size(orbs, 1)) error stop "SysXc_AddDensity: density and orbitals differ in size"

    !$omp parallel do default(shared) private(p, l, acc) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do p = 1, size(density)
      acc = 0.0_R64
      do l = 1, size(orbs, 2)
        acc = acc + real(orbs(p, l), kind=R64) * real(weightedOrbs(p, l), kind=R64) + &
              aimag(orbs(p, l)) * aimag(weightedOrbs(p, l))
      end do
      density(p) = density(p) + acc
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> resOrb_p = pot_p orb_p (diagonal metric: Hermitian for real pot).
  subroutine MultiplyWithPotentialNodal(resOrb, pot, orb)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(out), contiguous, target :: resOrb(:)
    real(R64), intent(in), contiguous, target :: pot(:)
    complex(R64), intent(in), contiguous, target :: orb(:)

    integer(I32) :: p

    !$omp parallel do default(shared) private(p) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do p = 1, size(resOrb)
      resOrb(p) = pot(p) * orb(p)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> src_p = ρ_p w_p (SysInteraction_FillInteractionSrc of the nodal families:
  !> conjg(a) b w).
  subroutine FillInteractionSrcNodal(src, density)
    complex(R64), intent(out), allocatable, target :: src(:)
    real(R64), intent(in), contiguous, target :: density(:)

    allocate (src(size(density)))
    src(:) = cmplx(density(:) * SysXc_weights(:), 0.0_R64, kind=R64)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> pot_p += Re V_p (nodal interaction potential).
  subroutine AddInteractionPotentialNodal(pot, interactionPotential)
    real(R64), intent(inout), contiguous, target :: pot(:)
    complex(R64), intent(in), contiguous, target :: interactionPotential(:)

    pot(:) = pot(:) + real(interactionPotential(:), kind=R64)

  end subroutine

end submodule
