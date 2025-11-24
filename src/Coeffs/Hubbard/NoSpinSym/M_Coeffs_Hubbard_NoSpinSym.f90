! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Hubbard backend: no spin symmetry (full direct-product space).
!>
!> Provides the most general coefficient representation as the direct product of
!> spin-up and spin-down configuration spaces. The total number of coefficients is
!> `n_up * n_dn`. This backend makes no use of (anti)symmetry under spin exchange
!> and is therefore the most flexible but also the largest state space.
module M_Coeffs_Hubbard_NoSpinSym
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the no-spin-symmetry Hubbard representation.
    !>
    !> - compute `Coeffs_nCoeffs = n_up * n_dn`,
    !> - define the (i_up, i_dn) <-> i linearization for the full product space,
    !> - attach implementations for H1/H2 application, RDM fills, index mapping,
    !>   normalization, projection, and persistence via `M_Coeffs` pointers.
    module subroutine Coeffs_Hubbard_NoSpinSym_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
