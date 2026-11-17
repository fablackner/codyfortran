! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Hubbard backend without spin symmetry: the full product space of the
!> configurations of the two spins, c(iUp, iDn) stored with the up index
!> fastest, n_up n_dn coefficients. Any particle numbers and any h1 (also
!> different for the two spins).
!> The Hamiltonian application, RDMs and index maps are those of
!> M_Coeffs_Hubbard, which reads the parity set here.
module M_Coeffs_Hubbard_NoSpinSym
  use M_Utils_Types
  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================
  interface
    !> Selects the product space: parity 0, n_up n_dn coefficients.
    module subroutine Coeffs_Hubbard_NoSpinSym_Fabricate
    end subroutine
  end interface

end module
