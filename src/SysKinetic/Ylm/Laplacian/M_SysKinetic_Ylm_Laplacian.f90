! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ylm radial Laplacian backends for kinetic operators.
!>
!> Hosts shared data and fabrication for ylm radial-kinetic implementations,
!> such as finite-difference and FEDVR. The operator targets individual (l,m)
!> channels and is used together with angular momentum terms.
module M_SysKinetic_Ylm_Laplacian
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate and configure the ylm radial Laplacian backend.
    !>
    !> Selects the numerical scheme (FinDiff/FEDVR) and prepares masses and
    !> any radial grid metadata needed by the chosen implementation.
    module subroutine SysKinetic_Ylm_Laplacian_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Body masses used for radial kinetic scaling via 1/(2 m_bt)
  real(R64), allocatable :: SysKinetic_Ylm_Laplacian_bodyMass(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
