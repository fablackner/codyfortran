! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Hubbard-model coefficient backend.
!>
!> This module provides data and runtime wiring for CI coefficients tailored to
!> the Hubbard model. It uses compact bit encodings for spin-resolved
!> configurations and precomputes connectivity and weights for efficient
!> application of the kinetic (hopping) and interaction terms.
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
