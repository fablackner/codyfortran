! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysKinetic_Ylm_Laplacian) S_SysKinetic_Ylm_Laplacian

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysKinetic_Ylm_Laplacian_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Ylm_Laplacian_FinDiff
    use M_SysKinetic_Ylm_Laplacian_Fedvr

    call Say_Fabricate("sysKinetic.ylm.laplacian")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysKinetic_Ylm_Laplacian_bodyMass = Json_Get("sysKinetic.ylm.laplacian.bodyMass", [1.0_R64])

    SysKinetic_timeIndependentQ = .true.
    SysKinetic_bodyTypeIndependentQ = .true.

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysKinetic.ylm.laplacian.finDiff")) then
      call SysKinetic_Ylm_Laplacian_FinDiff_Fabricate

    else if (Json_GetExistence("sysKinetic.ylm.laplacian.fedvr")) then
      call SysKinetic_Ylm_Laplacian_Fedvr_Fabricate

    else
      error stop "sysKinetic.ylm.laplacian is missing one of: finDiff, fedvr"
    end if
  end subroutine

end submodule
