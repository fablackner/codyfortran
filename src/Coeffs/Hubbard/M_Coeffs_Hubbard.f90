! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Hubbard-model coefficient backend with bit-encoded Fock states.
!>
!> # Purpose
!>
!> This module provides an optimized CI representation for the Fermi-Hubbard model
!> on lattice grids. Instead of the general tensor-product enumeration used by
!> Generic, configurations are encoded as bit patterns where bit j=1 indicates
!> site j is occupied.
!>
!> # Bit Encoding
!>
!> For N_orb lattice sites and N_up (N_dn) up-spin (down-spin) electrons:
!> - Each spin-up config is a 64-bit integer with exactly N_up bits set
!> - Total up-configs: C(N_orb, N_up) stored in `bitcodesUP(:)`
!> - Total down-configs: C(N_orb, N_dn) stored in `bitcodesDN(:)`
!>
!> # Hopping Graphs
!>
!> Single-electron hops (kinetic term) are precomputed during `Setup`:
!> - `hoppUP(k, i)`: Target config index when applying k-th hop to config i
!> - `weightUP(k, i)`: Matrix element including fermionic sign from bit counting
!> - `nConnectedUP(i)`: Number of valid hops from config i
!>
!> # Spin Symmetry Variants
!>
!> Three sub-backends exploit spin exchange symmetry:
!>
!> | Backend        | State space dimension | Physical interpretation          |
!> |----------------|----------------------|----------------------------------|
!> | `noSpinSym`    | n_up × n_dn          | Full product space               |
!> | `plusSpinSym`  | n(n+1)/2             | Symmetric (triplet-like) sector  |
!> | `minusSpinSym` | n(n-1)/2             | Antisymmetric (singlet-like)     |
!>
!> # JSON Configuration
!>
!> ```json
!> "coeffs": {
!>   "hubbard": {
!>     "noSpinSym": { }     // or "plusSpinSym": { } or "minusSpinSym": { }
!>   }
!> }
!> ```
!>
!> # Performance Notes
!>
!> - Bit operations (`popcnt`, `btest`, `ibset`, `ibclr`) enable O(1) sign computation
!> - Hopping graph is sparse: only ~2×nDim hops per config (nDim = lattice dimensions)
!> - Interaction is diagonal: on-site U·n↑n↓ via `interactionValues` lookup
!>
!> @see M_Coeffs_Generic     Alternative for non-Hubbard models
!> @see M_Grid_Lattice       Lattice geometry definitions
!> @see M_SysKinetic_Lattice Hopping matrix construction
module M_Coeffs_Hubbard
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Select and initialize the Hubbard coefficient representation.
    !>
    !> Determines array sizes and assigns `M_Coeffs` procedure pointers according to
    !> the chosen spin-symmetry backend (no spin symmetry, plus symmetry, minus
    !> symmetry). May delegate to the symmetry-specific modules.
    module subroutine Coeffs_Hubbard_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Number of distinct spin-up configurations.
  !> Equals the binomial coefficient C(n_orb, n_up).
  integer(I32) :: Coeffs_Hubbard_nCoeffsUP

  !> Number of distinct spin-down configurations.
  !> Equals the binomial coefficient C(n_orb, n_dn).
  integer(I32) :: Coeffs_Hubbard_nCoeffsDN

  !> For each up-configuration i, number of single-hop target configurations.
  !> Used to size slices in the `hoppUP`/`weightUP` tables.
  integer(I32), allocatable :: Coeffs_Hubbard_nConnectedUP(:)

  !> For each down-configuration i, number of single-hop target configurations.
  !> Used to size slices in the `hoppDN`/`weightDN` tables.
  integer(I32), allocatable :: Coeffs_Hubbard_nConnectedDN(:)

  !> Hopping graph (spin up): target indices for single-electron hops.
  !> `hoppUP(k,i)` gives the index (1..nCoeffsUP) of the k-th target reachable
  !> from source configuration `i` by one allowed hop. Only the first
  !> `nConnectedUP(i)` entries in column `i` are valid.
  integer(I32), allocatable :: Coeffs_Hubbard_hoppUP(:, :)

  !> Hopping graph (spin down): target indices for single-electron hops.
  !> `hoppDN(k,i)` gives the index (1..nCoeffsDN) of the k-th target reachable
  !> from source configuration `i` by one allowed hop. Only the first
  !> `nConnectedDN(i)` entries in column `i` are valid.
  integer(I32), allocatable :: Coeffs_Hubbard_hoppDN(:, :)

  !> Hopping amplitudes (spin up): matrix elements for the kinetic operator.
  !> `weightUP(k,i)` corresponds to the transition `i -> hoppUP(k,i)`.
  real(R64), allocatable :: Coeffs_Hubbard_weightUP(:, :)

  !> Hopping amplitudes (spin down): matrix elements for the kinetic operator.
  !> `weightDN(k,i)` corresponds to the transition `i -> hoppDN(k,i)`.
  real(R64), allocatable :: Coeffs_Hubbard_weightDN(:, :)

  !> Bit representation of spin-up configurations (length = n_orb bits used).
  !> Bit j set means lattice site j is occupied by an up electron.
  integer(I64), allocatable :: Coeffs_Hubbard_bitcodesUP(:)

  !> Bit representation of spin-down configurations (length = n_orb bits used).
  !> Bit j set means lattice site j is occupied by a down electron.
  integer(I64), allocatable :: Coeffs_Hubbard_bitcodesDN(:)

  !> Onsite interaction prefactors for each (up,down) pair in the direct-product
  !> space. Typically derived from popcount(bitcodesUP(i) AND bitcodesDN(j)).
  !> The linearization order is consistent with the representation’s
  !> index mapping used by `ConfigurationsFromIndex`.
  real(R64), allocatable :: Coeffs_Hubbard_interactionValues(:)

  !=============================================================================
  ! module procedures pointers
  ! (exported via `M_Coeffs` by the symmetry backends)
  !=============================================================================

end module
