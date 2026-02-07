! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Lattice_Manual) S_SysPotential_Lattice_Manual

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Lattice_Manual_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Lattice_Manual_StdImpl

    implicit none

    real(R64), allocatable :: values(:)
    integer(I32), allocatable :: sites(:)

    call Say_Fabricate("sysPotential.lattice.manual")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    if (Json_GetExistence("sysPotential.lattice.manual.values")) then
      values = Json_Get("sysPotential.lattice.manual.values", [0.0_R64])
      SysPotential_Lattice_Manual_values = values
    else
      error stop "sysPotential.lattice.manual is missing: values"
    end if

    if (Json_GetExistence("sysPotential.lattice.manual.sites")) then
      sites = Json_Get("sysPotential.lattice.manual.sites", [1])
      SysPotential_Lattice_Manual_sites = sites
    end if

    SysPotential_timeIndependentQ = .true.
    SysPotential_bodyTypeIndependentQ = .true.

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.lattice.manual.stdImpl")) then
      call SysPotential_Lattice_Manual_StdImpl_Fabricate

    else
      error stop "sysPotential.lattice.manual is missing one of: stdImpl"
    end if

  end subroutine

end submodule
