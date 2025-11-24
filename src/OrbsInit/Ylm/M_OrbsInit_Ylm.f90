! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Spherical-harmonics (Y_l^m) orbital initialization backend.
!>
!> Provides an initialization function interface for orbitals defined in a
!> spherical basis. Typical use is to construct radial-times-angular products
!> where the angular part is Y_l^m. The fabricate routine connects Ylm-specific
!> implementations to the generic pointers in `M_OrbsInit`.
module M_OrbsInit_Ylm
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Select and wire the Ylm-based initialization implementation.
    !> Reads any Ylm-related configuration (quantum numbers, normalization,
    !> radial model) and assigns the `M_OrbsInit` procedure pointers.
    module subroutine OrbsInit_Ylm_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointwise initializer for spherical-harmonics based orbitals.
  procedure(I_OrbsInit_Ylm_InitFunction), pointer :: OrbsInit_Ylm_InitFunction
  abstract interface
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !> Compute the complex orbital value at radius r for given (l, m, index).
    function I_OrbsInit_Ylm_InitFunction(r, l, m, index, bt_) result(res)
      import :: I32, R64
      !> Resulting complex orbital value.
      complex(R64) :: res
      !> Radial coordinate r (units consistent with backend).
      real(R64), intent(in) :: r
      !> Orbital angular momentum quantum number l.
      integer(I32), intent(in) :: l
      !> Magnetic quantum number m.
      integer(I32), intent(in) :: m
      !> Orbital index within the current species/body type.
      integer(I32), intent(in) :: index
      !> Optional body type/species/bath index.
      integer(I32), intent(in), optional :: bt_
    end function
  end interface

end module
