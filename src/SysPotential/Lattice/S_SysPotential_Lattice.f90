! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Lattice-grid fabrication submodule for external potentials.
!>
!> Dispatches to the appropriate lattice potential model (harmonic, disorder,
!> seesaw, manual) based on the JSON configuration. All lattice potentials share
!> the same point-wise multiplication operator implemented here.
submodule(M_SysPotential_Lattice) S_SysPotential_Lattice

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Parse lattice potential configuration and wire the selected model.
  !>
  !> Available models:
  !> - `harmonic`:      Parabolic trap V = ½ω²(r-r₀)²
  !> - `randomUniform`: Disorder from uniform distribution
  !> - `randomGauss`:   Disorder from Gaussian distribution
  !> - `seesaw`:        Time-dependent linear tilt
  !> - `manual`:        User-specified site values
  module subroutine SysPotential_Lattice_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Lattice_Harmonic
    use M_SysPotential_Lattice_RandomUniform
    use M_SysPotential_Lattice_RandomGauss
    use M_SysPotential_Lattice_Seesaw
    use M_SysPotential_Lattice_Manual

    call Say_Fabricate("sysPotential.lattice")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_MultiplyWithExternalPotential => MultiplyWithExternalPotential

    !------------------------------------
    ! branch by potential model
    !------------------------------------

    if (Json_GetExistence("sysPotential.lattice.harmonic")) then
      call SysPotential_Lattice_Harmonic_Fabricate

    else if (Json_GetExistence("sysPotential.lattice.randomUniform")) then
      call SysPotential_Lattice_RandomUniform_Fabricate

    else if (Json_GetExistence("sysPotential.lattice.randomGauss")) then
      call SysPotential_Lattice_RandomGauss_Fabricate

    else if (Json_GetExistence("sysPotential.lattice.seesaw")) then
      call SysPotential_Lattice_Seesaw_Fabricate

    else if (Json_GetExistence("sysPotential.lattice.manual")) then
      call SysPotential_Lattice_Manual_Fabricate

    else
      error stop "sysPotential.lattice is missing one of: harmonic, randomUniform, randomGauss, seesaw, manual"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Apply the external potential to an orbital via point-wise multiplication.
  !>
  !> For lattice grids, the potential is diagonal in the site basis, so
  !> application is simply: dOrb(i) = V(i) × orb(i).
  subroutine MultiplyWithExternalPotential(dOrb, externalPotential, orb)
    use M_Grid

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: externalPotential(:)
    complex(R64), intent(in), contiguous :: orb(:)

    dOrb(:) = externalPotential(:) * orb(:)

  end subroutine

end submodule
