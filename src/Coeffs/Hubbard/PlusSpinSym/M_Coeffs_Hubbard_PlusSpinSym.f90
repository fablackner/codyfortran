! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Hubbard backend in the sector symmetric under the exchange of the two
!> spins: the orthonormal basis (|i,k> + |k,i>)/sqrt(2) (i < k) and |i,i> of
!> the configurations i, k of one spin, n(n+1)/2 coefficients stored by
!> columns of the upper triangle (index (k - 1) k / 2 + i for i <= k). Needs
!> equal particle numbers of the two spins and a spin-independent h1.
!> The Hamiltonian application, RDMs and index maps are those of
!> M_Coeffs_Hubbard, which reads the parity set here.
module M_Coeffs_Hubbard_PlusSpinSym
  use M_Utils_Types
  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================
  interface
    !> Selects the symmetric sector: parity +1, n(n+1)/2 coefficients.
    module subroutine Coeffs_Hubbard_PlusSpinSym_Fabricate
    end subroutine
  end interface

end module
