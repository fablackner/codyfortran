! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Static nuclear Coulomb attraction of the grid.fem centers.
!>
!>   V(r) = sum_I V_I(r),   nucleus I at R_I = Grid_Fem_centers(:, I),
!>
!> with, per center (charge Z_I, smearing radius rc_I):
!>   - rc_I = 0: Galerkin point nucleus. All point nuclei together are ONE
!>     Grid_Fem_SolvePoissonPointCharges solve with the charges q_I = -Z_I:
!>     load_i = sum_I q_I phi_i(R_I), exact box values sum_I q_I/|r_B - R_I|.
!>     The nodal potential is -Z/|r - R| up to the mesh error away from the
!>     nucleus and finite at the nucleus node (roughly -30 ... -50 Z/hMin
!>     for nLocals 5-6, depending on the mesh, about -85 Z/hMin for
!>     nLocals 8; hMin the local element size).
!>     The electron-nucleus attraction sum_i w_i rho_i V_i then equals -Z
!>     times the Galerkin Hartree potential of rho at R_I, i.e. it is as
!>     accurate as the FEM Hartree potential.
!>   - rc_I > 0: smeared nucleus, the nodal values -Z erf(|r - R|/rc)/|r - R|
!>     (limit -2Z/(sqrt(pi) rc) at the nucleus node). Smooth, no spike: the
!>     physics changes in the core only, but the spectrum of T + V is no longer
!>     stretched by the nuclear spike, which allows explicit time steps of
!>     about 0.01-0.05 a.u. on correspondingly coarser core meshes.
!>   - Z_I = 0: ghost center (mesh refinement only, no potential).
!>
!> The potential is time- and body-type-independent; Setup builds it once
!> (a private cache of the implementation submodule) and FillExternalPotential
!> copies it. The charges, smearing radii and the nuclear repulsion
!> sum_(I<J) Z_I Z_J/R_IJ (printed at Setup) are module data, for programs
!> and printers.
!>
!> JSON Configuration
!> ------------------
!>   "sysPotential": { "fem": { "coulomb": {
!>       "charges":  [8.0, 1.0, 1.0],     // required, one per grid.fem center
!>       "smearing": [0.0, 0.0, 0.0]      // optional, default 0 (point nuclei)
!>   } } }
module M_SysPotential_Fem_Coulomb
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Read and validate charges and smearing radii, compute the nuclear
    !> repulsion, set the flags and bind Setup and FillExternalPotential.
    module subroutine SysPotential_Fem_Coulomb_Fabricate
    end subroutine

    !> Potential -Z erf(r/rc)/r of a smeared nucleus at distance r >= 0 (the
    !> limit -2Z/(sqrt(pi) rc) for r -> 0; no division by r at the nucleus).
    pure module function SysPotential_Fem_Coulomb_SmearedPotential(charge, smearing, r) result(res)
      !> Nuclear charge Z.
      real(R64), intent(in) :: charge
      !> Smearing radius rc > 0.
      real(R64), intent(in) :: smearing
      !> Distance r >= 0 from the nucleus.
      real(R64), intent(in) :: r
      real(R64) :: res
    end function
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Nuclear charge Z_I of every grid.fem center (nCenters entries; 0 marks a
  !> ghost center, which only refines the mesh).
  real(R64), allocatable :: SysPotential_Fem_Coulomb_charges(:)

  !> Smearing radius rc_I of every center (nCenters entries): 0 for a Galerkin
  !> point nucleus, rc > 0 for the smeared nucleus -Z erf(r/rc)/r.
  real(R64), allocatable :: SysPotential_Fem_Coulomb_smearing(:)

  !> Nuclear repulsion sum_(I<J) Z_I Z_J / R_IJ of the point charges (set at
  !> fabrication; not part of the electronic energies of Method).
  real(R64) :: SysPotential_Fem_Coulomb_nuclearRepulsion = 0.0_R64

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
