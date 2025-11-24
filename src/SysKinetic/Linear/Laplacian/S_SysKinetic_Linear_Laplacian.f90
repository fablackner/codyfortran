! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysKinetic_Linear_Laplacian) S_SysKinetic_Linear_Laplacian

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysKinetic_Linear_Laplacian_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Linear_Laplacian_FinDiff
    use M_SysKinetic_Linear_Laplacian_Fourier

    call Say_Fabricate("sysKinetic.linear.laplacian")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysKinetic_Linear_Laplacian_bodyMass = Json_Get("sysKinetic.linear.laplacian.bodyMass", [1.0_R64])

    SysKinetic_timeIndependentQ = .true.
    SysKinetic_bodyTypeIndependentQ = .true.

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysKinetic.linear.laplacian.finDiff")) then
      call SysKinetic_Linear_Laplacian_FinDiff_Fabricate

    else if (Json_GetExistence("sysKinetic.linear.laplacian.fourier")) then
      call SysKinetic_Linear_Laplacian_Fourier_Fabricate

    else
      error stop "sysKinetic.linear.laplacian is missing one of: finDiff"
    end if
  end subroutine

end submodule
