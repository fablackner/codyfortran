! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic.f90
!> @brief Implementation submodule for SysKinetic fabrication.
!>
!> @details
!> Dispatches to the appropriate backend based on JSON configuration:
!>
!>   - `sysKinetic.lattice` → M_SysKinetic_Lattice (tight-binding hopping)
!>   - `sysKinetic.linear`  → M_SysKinetic_Linear (Cartesian Laplacian)
!>   - `sysKinetic.ylm`     → M_SysKinetic_Ylm (spherical radial kinetic)
!>
!> The fabrication pattern ensures exactly one backend is active. Each backend
!> then binds the global `SysKinetic_MultiplyWithKineticOp` and `SysKinetic_Setup`
!> procedure pointers.
submodule(M_SysKinetic) S_Kinetic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Select and configure the kinetic backend from JSON.
  !>
  !> @details
  !> Examines the JSON tree under "sysKinetic" and delegates to the matching
  !> backend fabricator. Exactly one of {lattice, linear, ylm} must be present.
  module subroutine SysKinetic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic_Lattice
    use M_SysKinetic_Linear
    use M_SysKinetic_Ylm

    call Say_Fabricate("sysKinetic")

    !------------------------------------
    ! branch to the appropriate backend
    !------------------------------------

    if (Json_GetExistence("sysKinetic.lattice")) then
      call SysKinetic_Lattice_Fabricate

    else if (Json_GetExistence("sysKinetic.linear")) then
      call SysKinetic_Linear_Fabricate

    else if (Json_GetExistence("sysKinetic.ylm")) then
      call SysKinetic_Ylm_Fabricate

    else
      error stop "sysKinetic is missing one of: lattice, linear, ylm"
    end if

  end subroutine

end submodule
