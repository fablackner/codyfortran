! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Hubbard-model coefficient backend with bit-encoded Fock states.
!>
!> # Purpose
!>
!> CI for two fermion species (spin up and spin down, body types 1 and 2) in
!> the same orbitals, whose two-body Hamiltonian is an on-site interaction
!>   H2 = sum_s U_s n_{s,up} n_{s,dn}
!> in the orbital basis: the Fermi-Hubbard model with the lattice sites as
!> orbitals (orbsInit.lattice.onSite, static orbitals: method tdci). The CI
!> vector is the matrix c(iUp, iDn) of the configurations of the two spins;
!> H1 acts on its rows and columns separately and H2 is diagonal.
!>
!> # Bit encoding
!>
!> A configuration of one spin is a 64-bit integer, bit o - 1 set if orbital o
!> is occupied (at most 62 orbitals). The configurations of N particles are
!> numbered in increasing order of their bit pattern (colex order), so the
!> index of a pattern is one plus its rank in the combinatorial number system,
!> sum_t binom(position of the t-th set bit, t): no lookup table of size 2^n.
!> Fermionic signs follow the orbital order: an operator on orbital o picks up
!> (-1)^(occupied orbitals below o).
!>
!> # Single excitations
!>
!> `Coeffs_Setup` tabulates a+_c a_d |i> = sign |target> for every
!> configuration i of each spin and every pair (c, d) with d occupied (c = d
!> included). The tables do not depend on the Hamiltonian: H1 is built from
!> the h1_ of each call (time-dependent h1, mean-field and basis-change
!> operators included), and the same tables give the 1-RDM and 2-RDM.
!>
!> # Spin symmetry variants
!>
!> | Backend        | Coefficients | Stored vector                                |
!> |----------------|--------------|----------------------------------------------|
!> | `noSpinSym`    | n_up n_dn    | c(iUp, iDn), up index fastest                |
!> | `plusSpinSym`  | n(n+1)/2     | (|i,k> + |k,i>)/sqrt(2) for i < k, |i,i>     |
!> | `minusSpinSym` | n(n-1)/2     | (|i,k> - |k,i>)/sqrt(2) for i < k            |
!>
!> The sectors are the eigenspaces of the exchange of the two spins; they need
!> equal particle numbers and a spin-independent h1. The basis of a sector is
!> orthonormal, so norms and inner products of the stored vectors are those of
!> the states.
!>
!> # JSON Configuration
!>
!> ```json
!> "coeffs": { "hubbard": { "noSpinSym": { } } }
!> ```
!> (or "plusSpinSym": { } or "minusSpinSym": { }); no configList is needed.
!>
!> @see M_Coeffs_Generic     General CI of any body types and Hamiltonians
module M_Coeffs_Hubbard
  use M_Utils_Types
  implicit none

  !=============================================================================
  ! types
  !=============================================================================

  !> Single excitations of the configurations of one spin, by source
  !> configuration: the entries first(i) .. first(i + 1) - 1 hold
  !> a+_create(k) a_destroy(k) |i> = factor(k) |target(k)> (orbital indices
  !> of the spin, 1 .. nSites; factor = +-1).
  type :: T_Coeffs_Hubbard_Singles
    integer(I32), allocatable :: first(:)
    integer(I32), allocatable :: target(:)
    integer(I32), allocatable :: create(:)
    integer(I32), allocatable :: destroy(:)
    real(R64), allocatable :: factor(:)
  end type

  !=============================================================================
  ! module procedures
  !=============================================================================
  interface
    !> Checks the many-body setting (two fermion species in the same number
    !> of orbitals), counts the configurations of each spin, binds the
    !> backend procedures of M_Coeffs and selects the spin-symmetry variant.
    module subroutine Coeffs_Hubbard_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Number of orbitals of each spin (the lattice sites).
  integer(I32) :: Coeffs_Hubbard_nSites

  !> Number of configurations of the spin-up particles, binom(nSites, N_up).
  integer(I32) :: Coeffs_Hubbard_nCoeffsUP

  !> Number of configurations of the spin-down particles, binom(nSites, N_dn).
  integer(I32) :: Coeffs_Hubbard_nCoeffsDN

  !> Stored representation, set by the variant: 0 the product space
  !> (noSpinSym), +1 or -1 the sector of that parity under the exchange of
  !> the spins (plusSpinSym, minusSpinSym).
  integer(I32) :: Coeffs_Hubbard_parity = 0

  !> Bit patterns of the spin-up configurations, in increasing order.
  integer(I64), allocatable :: Coeffs_Hubbard_bitcodesUP(:)

  !> Bit patterns of the spin-down configurations, in increasing order.
  integer(I64), allocatable :: Coeffs_Hubbard_bitcodesDN(:)

  !> Single excitations of the spin-up configurations.
  type(T_Coeffs_Hubbard_Singles) :: Coeffs_Hubbard_singlesUP

  !> Single excitations of the spin-down configurations.
  type(T_Coeffs_Hubbard_Singles) :: Coeffs_Hubbard_singlesDN

end module
