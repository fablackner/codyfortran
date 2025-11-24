! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential) S_Potential

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_Characters
    use M_Utils_NoOpProcedures
    use M_SysPotential_Lattice
    use M_SysPotential_Linear
    use M_SysPotential_Ylm

    call Say_Fabricate("sysPotential")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.lattice")) then
      call SysPotential_Lattice_Fabricate

    else if (Json_GetExistence("sysPotential.linear")) then
      call SysPotential_Linear_Fabricate

    else if (Json_GetExistence("sysPotential.ylm")) then
      call SysPotential_Ylm_Fabricate

    else
      error stop "sysPotential is missing one of: lattice, linear, ylm"

    end if

  end subroutine

end submodule
