! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Orbital API and runtime dispatch for the orbital part of the wavefunction.
!>
!> M_Orbs exposes a small, stable front-end API while the actual
!> implementation is selected at runtime. All operational procedures are
!> procedure pointers that get bound by `Orbs_Fabricate`, allowing different
!> orbital representations/backends (e.g., restricted vs. unrestricted) to be
!> swapped without touching call sites.
!>
!> Responsibilities
!> - Hold global orbital state: `Orbs_orbs`, `Orbs_nOrbsInState`.
!> - Provide procedure-pointer hooks for setup, projection, orthonormalization,
!>   and persistence which are bound to a concrete backend at runtime.
!>
!> Typical lifecycle
!> 1. Call `Orbs_Fabricate()` once to select and bind the backend.
!> 2. Call `Orbs_Setup()` (bound by the backend) to allocate/associate
!>    `Orbs_orbs` and initialize its contents.
!> 3. Use `Orbs_ProjectOnSubspace`, `Orbs_Orthonormalize`, and `Orbs_SaveOrbs`
!>    as needed during propagation and I/O.
!>
!> Notes
!> - Thread-safety: this module maintains global state and is not inherently
!>   thread-safe.
!> - Error handling and persistence details are backend-specific; see the
!>   concrete implementation selected by `Orbs_Fabricate`.
module M_Orbs
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Select and bind the active orbital backend and configure global flags.
    !>
    !> This routine must be called exactly once before using any `Orbs_*`
    !> procedures. It records whether a spin-restricted representation is
    !> requested and binds all orbital-related procedure pointers to the
    !> selected backend implementation. Memory is not allocated here; call
    !> `Orbs_Setup()` afterwards to allocate/associate `Orbs_orbs` and finalize
    !> initialization.
    module subroutine Orbs_Fabricate()
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Number of active single-particle orbitals in the current representation.
  !>
  !> Semantics: equals the number of columns in `Orbs_orbs` once setup is
  !> complete. Backends are responsible for setting this value during
  !> `Orbs_Setup()`.
  integer(I32)                      :: Orbs_nOrbsInState

  !> Orbital coefficient matrix (associated by the backend during setup).
  !>
  !> Shape and layout: `Orbs_orbs(nBasis, nOrbs)` where
  !> - `nBasis` is the dimension of the one-particle basis used by the backend,
  !> - `nOrbs`  equals `Orbs_nOrbsInState`.
  !>
  !> Convention: `Orbs_orbs(i, j)` is the i-th basis component of the j-th
  !> orbital. Columns are typically assumed to be orthonormal after
  !> `Orbs_Orthonormalize` or backend initialization.
  complex(R64), contiguous, pointer :: Orbs_orbs(:, :)

  !> If true, a spin/body-type restricted representation is used and a single
  !> set of orbitals is propagated for all body types. Otherwise the backend
  !> may differentiate orbitals per body type.
  logical                           :: Orbs_restrictedQ

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Backend-provided setup routine.
  !>
  !> Responsibilities (backend contract)
  !> - allocate/associate `Orbs_orbs` and set `Orbs_nOrbsInState`.
  !> - initialize orbitals (e.g., from guess, file, or analytic forms).
  !> - ensure consistency with `Orbs_restrictedQ`.
  !> - leave columns reasonably orthonormal or be prepared for a subsequent
  !>   call to `Orbs_Orthonormalize`.
  procedure(I_Orbs_Setup), pointer :: Orbs_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize and allocate/associate the global orbital storage.
    !>
    !> Side effects
    !> - Associates `Orbs_orbs` and sets `Orbs_nOrbsInState`.
    !> - May read input and/or write diagnostic information depending on the
    !>   backend configuration.
    subroutine I_Orbs_Setup
    end subroutine
  end interface

  !> Project a set of vectors onto the orthogonal complement of a given
  !> orbital subspace.
  !>
  !> Use this to enforce gauge conditions or to remove components along a
  !> reference set of orbitals. The reference orbitals are typically assumed to
  !> have orthonormal columns.
  procedure(I_Orbs_ProjectOnSubspace), pointer :: Orbs_ProjectOnSubspace
  abstract interface
    !> Orthogonal projection onto the complement of span(orbs).
    !>
    !> Math
    !> Given vectors (columns) in `dOrbs` and an orthonormal set `orbs`, this
    !> routine replaces each column `ψ` in `dOrbs` with
    !> `ψ' = ψ - sum_j <φ_j|ψ> φ_j`, i.e., removes components parallel to the
    !> subspace spanned by `orbs`.
    !>
    !> Requirements
    !> - `orbs` columns should be orthonormal or well-conditioned.
    !> - `dOrbs` and `orbs` must share the same leading dimension (basis size).
    subroutine I_Orbs_ProjectOnSubspace(dOrbs, orbs)
      import :: R64
      !> Input/output columns to be projected. On exit, each column holds the
      !> component orthogonal to the linear hull of `orbs`.
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      !> Reference orbitals defining the subspace to project out.
      complex(R64), intent(in), contiguous    :: orbs(:, :)
    end subroutine
  end interface

  !> Orthonormalize the columns of an orbital coefficient matrix.
  !>
  !> Implementations typically use modified Gram–Schmidt with (optional)
  !> re-orthogonalization or an equivalent stable method. The routine updates
  !> the input in place and is expected to produce
  !> <φ_i | φ_j> = δ_ij within numerical tolerance.
  procedure(I_Orbs_Orthonormalize), pointer :: Orbs_Orthonormalize
  abstract interface
    !> In-place orthonormalization of orbital columns.
    !>
    !> On exit, columns are mutually orthonormal with respect to the standard
    !> complex inner product on the one-particle basis.
    subroutine I_Orbs_Orthonormalize(orbs)
      import :: R64
      !> Input/output orbital matrix `orbs(nBasis, nOrbs)`.
      !> On exit, columns are orthonormalized in place.
      complex(R64), intent(inout), contiguous :: orbs(:, :)
    end subroutine
  end interface

  !> Persist the current set of orbitals to disk for analysis, restart, or I/O
  !> with external tools. The exact file format, naming, and metadata are
  !> backend-defined (e.g., one file per orbital or a single container file).
  procedure(I_Orbs_SaveOrbs), pointer :: Orbs_SaveOrbs
  abstract interface
    !> Save an orbital coefficient matrix to one or more files.
    !>
    !> The routine does not modify its input. Paths, formats, and auxiliary
    !> metadata are controlled by the backend and project-wide settings.
    subroutine I_Orbs_SaveOrbs(orbs)
      import :: R64
      !> Orbital matrix to save. Shape is implementation-dependent but
      !> typically `orbs(nBasis, nOrbs)`.
      complex(R64), intent(in), contiguous :: orbs(:, :)
    end subroutine
  end interface

end module
