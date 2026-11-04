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
!>   - `sysKinetic.prolate` → M_SysKinetic_Prolate (two-center Laplacian)
!>   - `sysKinetic.fem`     → M_SysKinetic_Fem (3D spectral-element Laplacian)
!>
!> The fabrication pattern ensures exactly one backend is active. Each backend
!> then binds the global `SysKinetic_MultiplyWithKineticOp` and `SysKinetic_Setup`
!> procedure pointers; the optional hook `SysKinetic_SolveShifted` is nullified
!> first and bound only by backends that provide it (Fem).
submodule(M_SysKinetic) S_Kinetic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Select and configure the kinetic backend from JSON.
  !>
  !> @details
  !> Examines the JSON tree under "sysKinetic" and delegates to the matching
  !> backend fabricator. Exactly one of {lattice, linear, ylm, prolate, fem}
  !> must be present.
  module subroutine SysKinetic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic_Lattice
    use M_SysKinetic_Linear
    use M_SysKinetic_Ylm
    use M_SysKinetic_Prolate
    use M_SysKinetic_Fem

    call Say_Fabricate("sysKinetic")

    ! Optional hooks are bound only by backends that provide them; reset them
    ! so that a backend fabricated earlier in the same process leaves no stale
    ! binding (e.g. the Fem shifted solve for a later non-Fem backend)
    nullify (SysKinetic_SolveShifted)

    !------------------------------------
    ! branch to the appropriate backend
    !------------------------------------

    if (Json_GetExistence("sysKinetic.lattice")) then
      call SysKinetic_Lattice_Fabricate

    else if (Json_GetExistence("sysKinetic.linear")) then
      call SysKinetic_Linear_Fabricate

    else if (Json_GetExistence("sysKinetic.ylm")) then
      call SysKinetic_Ylm_Fabricate

    else if (Json_GetExistence("sysKinetic.prolate")) then
      call SysKinetic_Prolate_Fabricate

    else if (Json_GetExistence("sysKinetic.fem")) then
      call SysKinetic_Fem_Fabricate

    else
      error stop "sysKinetic is missing one of: lattice, linear, ylm, prolate, fem"
    end if

  end subroutine

end submodule
