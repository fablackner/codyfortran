! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Fem.f90
!> @brief Implementation submodule for FEM kinetic fabrication.
submodule(M_SysKinetic_Fem) S_SysKinetic_Fem

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Validate the grid and dispatch to the kinetic model.
  module subroutine SysKinetic_Fem_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic_Fem_Laplacian

    call Say_Fabricate("sysKinetic.fem")

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.fem")) then
      error stop "grid.fem is required for sysKinetic.fem"
    end if

    !------------------------------------
    ! dispatch to the kinetic model
    !------------------------------------

    if (Json_GetExistence("sysKinetic.fem.laplacian")) then
      call SysKinetic_Fem_Laplacian_Fabricate

    else
      error stop "sysKinetic.fem is missing one of: laplacian"
    end if

  end subroutine

end submodule
