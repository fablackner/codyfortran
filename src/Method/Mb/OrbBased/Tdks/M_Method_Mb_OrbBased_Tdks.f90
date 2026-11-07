! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Time-dependent Kohn-Sham method (real-time TDDFT, `method.mb.orbBased.tdks`).
!>
!> The state is a set of Kohn-Sham orbitals (one Slater determinant of
!> fermions, occupations from the static CI vector of `coeffsInit.unary`, as
!> for TDHX); each orbital φ_i of spin s evolves with
!>
!>   i dφ_i/dt = (h + v_KS,s[ρ↑, ρ↓] - a K_s) φ_i,
!>
!> h = T + V_ext (+ the SysGauge laser), v_KS,s = J + v_xc,s the local
!> Kohn-Sham potential of SysXc (adiabatic functional of the instantaneous
!> spin densities) and a K_s the fraction a = SysXc_exactExchangeFraction of
!> exact exchange of a hybrid within spin s. The Kohn-Sham energy
!>
!>   E = Σ_i ⟨g_i|h φ_i⟩ + E_H + E_xc - a E_x,   E_x = ½ Σ_s Σ_i ⟨g_i|K_s φ_i⟩,
!>
!> (g = φ R, R the 1-RDM blocks) is `Method_GetEnergy`; the orbital equations
!> are its Hamiltonian flow, so it is conserved without a laser.
!>
!> Spin channels: restricted orbitals (`orbs.restrictedQ`) must be closed
!> shell (equal 1-RDM spin blocks) and use the total density (unpolarized
!> functional); unrestricted orbitals use the polarized functional with
!> ρ↑, ρ↓ of body types 1 and 2 (one body type: ρ↓ = 0).
!>
!> Split-step propagation as MCTDHX: `Method_Mb_OrbBased_TimeDerivativeOrbsLin`
!> (-i h φ, or only the gauge coupling with `exactFieldFreeQ`) and
!> `Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin` (-i (v_KS - a K) φ);
!> with `exactFieldFreeQ` the exact part (Ylm: eigen-decomposition with phi
!> functions; grid.fem: Cayley step) propagates h0 = T + V_ext. Every orbital
!> derivative is projected below `exactFieldFreeCut` (Ylm; no-op without
!> cut): a GGA split step is stable only if every propagated h0 mode has
!> E dt < pi (A_Method.md, "Time step with a GGA").
module M_Method_Mb_OrbBased_Tdks
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds the TDKS procedure pointers (setup, time derivatives, energy, exact part).
    module subroutine Method_Mb_OrbBased_Tdks_Fabricate
    end subroutine

    !> dOrbs += (v_KS,s - a K_s) orbs for the stored orbitals of the 1-RDM rdm1
    !> (nOrbsSum layout), without the one-body operator; the Kohn-Sham energies
    !> of these orbitals are stored in the module variables below.
    module subroutine Method_Mb_OrbBased_Tdks_ApplyKohnShamOp(dOrbs, orbs, rdm1, time)
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      complex(R64), intent(in), contiguous :: orbs(:, :)
      complex(R64), intent(in), contiguous :: rdm1(:, :)
      real(R64), intent(in) :: time
    end subroutine

    !> Local Kohn-Sham potentials v_KS,s = J + v_xc,s of the stored orbitals
    !> and 1-RDM on the xc points of SysXc, (SysXc_nPoints, nSpin), and the
    !> column s of every body type (spinOfBodyType); the local part of the
    !> frozen mean-field reference (Method_Mb_OrbBased_BuildFieldFreeMeanField).
    module subroutine Method_Mb_OrbBased_Tdks_FillKohnShamPotential(kohnShamPot, spinOfBodyType, orbs, rdm1, time)
      real(R64), intent(out), allocatable :: kohnShamPot(:, :)
      integer(I32), intent(out), allocatable :: spinOfBodyType(:)
      complex(R64), intent(in), contiguous :: orbs(:, :)
      complex(R64), intent(in), contiguous :: rdm1(:, :)
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Spin channels of the functional: 1 (restricted, total density) or 2
  integer(I32) :: Method_Mb_OrbBased_Tdks_nSpin = 0

  !> Energies of the last Kohn-Sham evaluation (Method_GetEnergy or a time
  !> derivative): Σ⟨g|hφ⟩ (Method_GetEnergy only), Hartree, exchange-
  !> correlation and the exact exchange E_x (enters as -a E_x)
  real(R64) :: Method_Mb_OrbBased_Tdks_oneBodyEnergy = 0.0_R64
  real(R64) :: Method_Mb_OrbBased_Tdks_hartreeEnergy = 0.0_R64
  real(R64) :: Method_Mb_OrbBased_Tdks_xcEnergy = 0.0_R64
  real(R64) :: Method_Mb_OrbBased_Tdks_exactExchangeEnergy = 0.0_R64

end module
