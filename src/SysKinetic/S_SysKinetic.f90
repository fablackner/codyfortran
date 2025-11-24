! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysKinetic) S_Kinetic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysKinetic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic_Lattice
    use M_SysKinetic_Linear
    use M_SysKinetic_Ylm

    call Say_Fabricate("sysKinetic")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysKinetic.lattice")) then
      call SysKinetic_Lattice_Fabricate

    else if (Json_GetExistence("sysKinetic.linear")) then
      call SysKinetic_Linear_Fabricate

    else if (Json_GetExistence("sysKinetic.ylm")) then
      call SysKinetic_Ylm_Fabricate

    else
      error stop "sysKinetic is missing one of: linear, ylm"
    end if

  end subroutine

end submodule
