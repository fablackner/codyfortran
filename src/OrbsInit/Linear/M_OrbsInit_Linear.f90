! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Linear-grid orbital initialization backend.
!>
!> Provides an initialization function interface for orbitals defined on a
!> 1D spatial coordinate x. The fabricate routine connects linear-specific
!> implementations (e.g., plane waves, harmonic eigenstates) to the generic
!> pointers in `M_OrbsInit`.
module M_OrbsInit_Linear
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Select and wire the linear-grid initialization implementation.
    !> Reads any linear-grid configuration (domain, spacing, basis parameters)
    !> and assigns the `M_OrbsInit` procedure pointers accordingly.
    module subroutine OrbsInit_Linear_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointwise initializer f(x; index, bt_) used by linear backends.
  procedure(I_OrbsInit_Linear_InitFunction), pointer :: OrbsInit_Linear_InitFunction
  abstract interface
    !> Compute the real-valued amplitude for a linear orbital at position x.
    function I_OrbsInit_Linear_InitFunction(x, index, bt_) result(res)
      import :: I32, R64
      !> Resulting orbital amplitude/value at x.
      real(R64)                             :: res
      !> Spatial coordinate.
      real(R64), intent(in)                 :: x
      !> Orbital index within the current species/body type.
      integer(I32), intent(in)              :: index
      !> Optional body type/species/bath index.
      integer(I32), intent(in), optional    :: bt_
    end function
  end interface

end module
