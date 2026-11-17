! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Hubbard backend in the sector antisymmetric under the exchange of the two
!> spins: the orthonormal basis (|i,k> - |k,i>)/sqrt(2) (i < k) of the
!> configurations i, k of one spin, n(n-1)/2 coefficients stored by columns
!> of the strict upper triangle (index (k - 2)(k - 1) / 2 + i for i < k).
!> Needs equal particle numbers of the two spins and a spin-independent h1.
!> The Hamiltonian application, RDMs and index maps are those of
!> M_Coeffs_Hubbard, which reads the parity set here.
module M_Coeffs_Hubbard_MinusSpinSym
  use M_Utils_Types
  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================
  interface
    !> Selects the antisymmetric sector: parity -1, n(n-1)/2 coefficients.
    module subroutine Coeffs_Hubbard_MinusSpinSym_Fabricate
    end subroutine
  end interface

end module
