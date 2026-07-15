! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Two-center nuclear Coulomb attraction for the prolate grid.
!>
!> Homonuclear diatomic with charge Z on each nucleus at the foci:
!>
!>   V(r) = -Z (1/r1 + 1/r2) = -(2 Z / a) xi / (xi^2 - eta^2)
!>
!> since r1 = a(xi+eta) and r2 = a(xi-eta). The (xi^2 - eta^2) denominator is
!> cancelled by the volume element in all matrix elements, so the nuclear
!> singularities are analytically regularized. The potential is
!> phi-independent (a pure m = 0 channel) and time-independent; the fixed
!> nuclear repulsion Z^2/R is a constant and is NOT included here.
!>
!> JSON Configuration
!> ------------------
!>   "sysPotential": { "prolate": { "coulomb": { "charge": 1.0 } } }
module M_SysPotential_Prolate_Coulomb
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the two-center Coulomb attraction.
    module subroutine SysPotential_Prolate_Coulomb_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Nuclear charge Z of each of the two (identical) nuclei.
  real(R64) :: SysPotential_Prolate_Coulomb_charge

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
