! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Linear-grid fabrication submodule for external potentials.
!>
!> Dispatches to the appropriate linear (1D real-space) potential model based on
!> the JSON configuration. All linear potentials share the same point-wise
!> multiplication operator implemented here.
submodule(M_SysPotential_Linear) S_SysPotential_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Parse linear potential configuration and wire the selected model.
  !>
  !> Available models:
  !> - `harmonic`:   Parabolic trap V = ½ω²(x-x₀)²
  !> - `softYukawa`: Multi-center softened Yukawa/Coulomb potential
  module subroutine SysPotential_Linear_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Linear_Harmonic
    use M_SysPotential_Linear_SoftYukawa

    call Say_Fabricate("sysPotential.linear")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_MultiplyWithExternalPotential => MultiplyWithExternalPotential

    !------------------------------------
    ! branch by potential model
    !------------------------------------

    if (Json_GetExistence("sysPotential.linear.harmonic")) then
      call SysPotential_Linear_Harmonic_Fabricate

    else if (Json_GetExistence("sysPotential.linear.softYukawa")) then
      call SysPotential_Linear_SoftYukawa_Fabricate

    else
      error stop "sysPotential.linear is missing one of: harmonic, softYukawa"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Apply the external potential to an orbital via point-wise multiplication.
  !>
  !> For linear grids, the potential is diagonal in position space, so
  !> application is simply: dOrb(i) = V(i) × orb(i).
  subroutine MultiplyWithExternalPotential(dOrb, externalPotential, orb)
    use M_Grid

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: externalPotential(:)
    complex(R64), intent(in), contiguous :: orb(:)

    dOrb(:) = externalPotential(:) * orb(:)

  end subroutine

end submodule
