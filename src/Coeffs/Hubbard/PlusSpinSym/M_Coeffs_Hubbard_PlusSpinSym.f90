! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Hubbard backend: plus spin symmetry (symmetric under spin exchange).
!>
!> Provides a reduced CI space for triplet-like sectors where exchanging spin-up
!> and spin-down configurations leaves the state invariant. If the underlying
!> per-spin configuration count is `n`, the number of unique coefficients becomes
!> \( n(n+1)/2 \). This module assigns the corresponding concrete procedures to
!> the pointers declared in `M_Coeffs` and defines the index mappings between
!> pair states and the reduced linear index.
module M_Coeffs_Hubbard_PlusSpinSym
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the plus-spin-symmetric Hubbard representation.
    !>
    !> - compute `Coeffs_nCoeffs = n(n+1)/2` from per-spin configuration counts,
    !> - define the (i_up, i_dn) -> i linearization consistent with symmetry,
    !> - attach implementations for H1/H2 application, RDM fills, index mapping,
    !>   normalization, projection, and persistence via `M_Coeffs` pointers.
    module subroutine Coeffs_Hubbard_PlusSpinSym_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
