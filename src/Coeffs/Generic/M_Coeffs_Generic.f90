! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Generic CI coefficient backend (model-agnostic tensor-product basis).
!>
!> # Purpose
!>
!> This module provides the coefficient representation for arbitrary many-body
!> systems that do not require specialized encodings. The CI space is constructed
!> as a tensor product of per-body-type configuration lists (from M_ConfigList):
!> \[
!>   |I\rangle = |C_1\rangle \otimes |C_2\rangle \otimes \cdots \otimes |C_{n_{bt}}\rangle
!> \]
!>
!> The total dimension is \(\prod_{bt} n_{configurations}(bt)\).
!>
!> # Features
!>
!> - **Flexible statistics:** Fermions, bosons, or mixtures via ConfigList
!> - **Efficient Hamiltonian application:** Uses precomputed `singles`/`doubles`
!>   connectivity graphs from ConfigList for sparse H1/H2 action
!> - **OpenMP parallelization:** Loop over CI indices with thread-local RDM accumulators
!> - **Linear addressing:** Multi-base encoding maps configuration tuples to linear indices
!>
!> # Index Convention
!>
!> A configuration tuple `(C_1, C_2, ..., C_nbt)` maps to linear index via:
!>   iCoeff = 1 + Σ_bt (C_bt - 1) × stride_bt
!> where stride_bt = Π_{bt' < bt} nConfigurations(bt').
!>
!> # JSON Configuration
!>
!> ```json
!> "coeffs": {
!>   "generic": { }
!> }
!> ```
!>
!> No additional parameters—the basis is determined entirely by `configList`.
!>
!> @see M_ConfigList     Configuration enumeration (Fock states)
!> @see M_Coeffs_Hubbard Alternative backend for Hubbard lattice models
module M_Coeffs_Generic
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Attach generic implementations to `M_Coeffs` and size the coefficient space.
    !>
    !> Responsibilities:
    !> - determine `Coeffs_nCoeffs` from the active configuration description,
    !> - associate procedure pointers in `M_Coeffs` with generic implementations,
    !> - perform any light-weight setup that is independent of a specific model.
    !>
    !> Heavy-weight, model-specific data (e.g., integrals) are not owned here and
    !> must be provided via arguments of the abstract interfaces when needed.
    module subroutine Coeffs_Generic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  ! (all exported via `M_Coeffs` after fabrication)
  !=============================================================================

end module

