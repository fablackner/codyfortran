! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Finite-difference radial Laplacian for ylm channels.
!>
!> Discretizes the radial part of the Laplacian on a (possibly nonuniform)
!> radial grid using finite differences. Special care is taken near r=0 to
!> respect regularity conditions and the l(l+1)/r^2 term.
module M_SysKinetic_Ylm_Laplacian_FinDiff
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Configure the finite-difference radial Laplacian and bind the operator.
    !>
    !> Expected inputs include the radial grid definition and FD order. The
    !> resulting operator acts per (l,m) channel and supports mass scaling.
    module subroutine SysKinetic_Ylm_Laplacian_FinDiff_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
