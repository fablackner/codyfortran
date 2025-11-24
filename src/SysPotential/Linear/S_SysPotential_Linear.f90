! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Linear) S_SysPotential_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    ! branch
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
  subroutine MultiplyWithExternalPotential(dOrb, externalPotential, orb)
    use M_Grid

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: externalPotential(:)
    complex(R64), intent(in), contiguous :: orb(:)

    dOrb(:) = externalPotential(:) * orb(:)

  end subroutine

end submodule
