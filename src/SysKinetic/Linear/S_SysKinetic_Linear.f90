! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Linear.f90
!> @brief Implementation submodule for linear-grid kinetic fabrication.
submodule(M_SysKinetic_Linear) S_SysKinetic_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Dispatch to the appropriate linear-grid scheme.
  !>
  !> Currently supports: `laplacian` (with FinDiff or Fourier sub-backends).
  module subroutine SysKinetic_Linear_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic_Linear_Laplacian

    call Say_Fabricate("sysKinetic.linear")

    !------------------------------------
    ! branch to the appropriate scheme
    !------------------------------------

    if (Json_GetExistence("sysKinetic.linear.laplacian")) then
      call SysKinetic_Linear_Laplacian_Fabricate

    else
      error stop "sysKinetic.linear is missing one of: laplacian"
    end if

  end subroutine

end submodule
