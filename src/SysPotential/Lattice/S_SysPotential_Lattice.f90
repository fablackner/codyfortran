! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Lattice) S_SysPotential_Lattice

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Lattice_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Lattice_Harmonic
    use M_SysPotential_Lattice_RandomUniform
    use M_SysPotential_Lattice_RandomGauss
    use M_SysPotential_Lattice_Seesaw

    call Say_Fabricate("sysPotential.lattice")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_MultiplyWithExternalPotential => MultiplyWithExternalPotential

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.lattice.harmonic")) then
      call SysPotential_Lattice_Harmonic_Fabricate

    else if (Json_GetExistence("sysPotential.lattice.randomUniform")) then
      call SysPotential_Lattice_RandomUniform_Fabricate

    else if (Json_GetExistence("sysPotential.lattice.randomGauss")) then
      call SysPotential_Lattice_RandomGauss_Fabricate

    else if (Json_GetExistence("sysPotential.lattice.seesaw")) then
      call SysPotential_Lattice_Seesaw_Fabricate

    else
      error stop "sysPotential.lattice is missing one of: harmonic, randomUniform, randomGauss, seesaw"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MultiplyWithExternalPotential(dOrb, externalPotential, orb)
    use M_Grid

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: externalPotential(:)
    complex(R64), intent(in), contiguous :: orb(:)

    dOrb(:) = externalPotential(:) * orb(:)

  end subroutine

end submodule
