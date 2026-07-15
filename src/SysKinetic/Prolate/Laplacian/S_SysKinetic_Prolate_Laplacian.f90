! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Prolate_Laplacian.f90
!> @brief Implementation submodule for prolate Laplacian fabrication.
submodule(M_SysKinetic_Prolate_Laplacian) S_SysKinetic_Prolate_Laplacian

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Read mass array and dispatch to the FEDVR backend.
  module subroutine SysKinetic_Prolate_Laplacian_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Prolate_Laplacian_Fedvr

    call Say_Fabricate("sysKinetic.prolate.laplacian")

    !------------------------------------
    ! read shared configuration
    !------------------------------------

    SysKinetic_Prolate_Laplacian_bodyMass = Json_Get("sysKinetic.prolate.laplacian.bodyMass", [1.0_R64])

    SysKinetic_timeIndependentQ = .true.
    SysKinetic_bodyTypeIndependentQ = .true.

    !------------------------------------
    ! dispatch to numerical method
    !------------------------------------

    if (Json_GetExistence("sysKinetic.prolate.laplacian.fedvr")) then
      call SysKinetic_Prolate_Laplacian_Fedvr_Fabricate

    else
      error stop "sysKinetic.prolate.laplacian is missing one of: fedvr"
    end if
  end subroutine

end submodule
