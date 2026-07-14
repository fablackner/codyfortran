! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ylm-grid fabrication submodule for external potentials.
!>
!> Dispatches to the appropriate spherical-harmonics potential model based on
!> the JSON configuration. Provides the generic `FillExternalPotential` that
!> assembles the full potential from (l,m) components, and `MultiplyWithExternalPotential`
!> that uses `Grid_Ylm_SpatialProduct` for proper angular coupling.
submodule(M_SysPotential_Ylm) S_SysPotential_Ylm

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Parse Ylm potential configuration and wire the selected model.
  !>
  !> Available models:
  !> - `coulomb`: Central Coulomb potential V(r) = -Z/r
  module subroutine SysPotential_Ylm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_SphericalHarmonics
    use M_Grid_Ylm
    use M_SysPotential
    use M_SysPotential_Ylm_Coulomb

    call Say_Fabricate("sysPotential.ylm")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_MultiplyWithExternalPotential => MultiplyWithExternalPotential
    SysPotential_FillExternalPotential => FillExternalPotential

    !------------------------------------
    ! branch by potential model
    !------------------------------------

    if (Json_GetExistence("sysPotential.ylm.coulomb")) then
      call SysPotential_Ylm_Coulomb_Fabricate

      ! Precompute Gaunt coefficients for the potential multiply
      ! (l1 <= lmaxPot, l2,l3 <= lmax)
      call SphericalHarmonics_EnsureGauntTable(max(Grid_Ylm_lmax, SysPotential_Ylm_lmax), &
                                               Grid_Ylm_lmax, &
                                               Grid_Ylm_lmax)

    else
      error stop "sysPotential.ylm is missing one of: coulomb"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Assemble the full external potential from its (l,m) components.
  !>
  !> Iterates over all (l,m) pairs up to `SysPotential_Ylm_lmax` and calls
  !> `SysPotential_Ylm_FillExternalPotentialRadial` for each component, then
  !> packs them into the flattened potential array using `Grid_Ylm_SetLmComponent`.
  subroutine FillExternalPotential(externalPotential, time, bt_)
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), allocatable :: externalPotential(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: l, m, potSize, nRad
    integer(I32) :: lmaxPot
    complex(R64), allocatable :: potLm(:)

    nRad = Grid_Ylm_nRadial
    lmaxPot = SysPotential_Ylm_lmax
    potSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    if (.not. allocated(externalPotential)) allocate (externalPotential(potSize))
    externalPotential = 0.0_R64

    allocate (potLm(nRad))

    do l = 0, lmaxPot
      do m = -l, l
        call SysPotential_Ylm_FillExternalPotentialRadial(potLm, l, m, time, bt_)
        call Grid_Ylm_SetLmComponent(externalPotential, l, m, potLm)
      end do
    end do

    deallocate (potLm)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Apply the external potential to an orbital via spherical-harmonics product.
  !>
  !> For Ylm grids, the product of two functions expanded in spherical harmonics
  !> requires Gaunt coefficients. This is handled by `Grid_Ylm_SpatialProduct`:
  !>   (V * ψ)_{lm}(r) = Σ_{l'm'; l''m''} G(...) V_{l'm'}(r) ψ_{l''m''}(r)
  subroutine MultiplyWithExternalPotential(dOrb, externalPotential, orb)
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: externalPotential(:)
    complex(R64), intent(in), contiguous :: orb(:)

    integer(I32) :: lmax, lmaxPot

    lmax = Grid_Ylm_lmax
    lmaxPot = SysPotential_Ylm_lmax

    call Grid_Ylm_SpatialProduct(dOrb, externalPotential, orb, lmax, lmaxPot, lmax)

  end subroutine

end submodule
