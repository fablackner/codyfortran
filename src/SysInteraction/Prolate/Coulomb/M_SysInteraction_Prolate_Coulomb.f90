! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Coulomb interaction for prolate grids.
!>
!> @details Pure 1/|r1 - r2| interaction, decomposed per azimuthal channel m
!> by expanding the eta dependence in normalized associated Legendre functions
!> (Neumann expansion). The expansion is truncated at degree
!> `sysInteraction.prolate.coulomb.lmax` (default and maximum nEta - 1, the
!> highest degree the eta quadrature resolves).
module M_SysInteraction_Prolate_Coulomb
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the Coulomb interaction for prolate grids.
    module subroutine SysInteraction_Prolate_Coulomb_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Coupling strength multiplying 1/|r1 - r2| (1.0 for electrons).
  real(R64) :: SysInteraction_Prolate_Coulomb_strength

  !> Truncation degree of the Legendre (Neumann) expansion in eta.
  integer(I32) :: SysInteraction_Prolate_Coulomb_lmax

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
