! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysKinetic_Linear) S_SysKinetic_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysKinetic_Linear_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic_Linear_Laplacian

    call Say_Fabricate("sysKinetic.linear")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysKinetic.linear.laplacian")) then
      call SysKinetic_Linear_Laplacian_Fabricate

    else
      error stop "sysKinetic.linear is missing one of: laplacian"
    end if

  end subroutine

end submodule
