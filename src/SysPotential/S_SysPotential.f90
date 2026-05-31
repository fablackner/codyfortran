! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Top-level fabrication submodule for external potentials.
!>
!> This submodule reads the JSON configuration to determine which grid family
!> (lattice, linear, or ylm) is requested and dispatches to the appropriate
!> family-level fabricator. The family fabricator in turn selects the concrete
!> potential model (harmonic, Coulomb, disorder, etc.) and wires the procedure
!> pointers exported by M_SysPotential.
submodule(M_SysPotential) S_Potential

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Parse the JSON configuration and dispatch to the appropriate grid backend.
  !>
  !> The JSON structure is:
  !> ```json
  !> { "sysPotential": { "<grid>": { "<model>": { ... } } } }
  !> ```
  !> where `<grid>` is one of: `lattice`, `linear`, `ylm`.
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
    ! branch by grid family
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
