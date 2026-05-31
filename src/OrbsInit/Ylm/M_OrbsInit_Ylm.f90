! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Spherical-harmonics (Y_lm) orbital initialization backend.
!>
!> @details
!> Provides an initialization interface for orbitals represented in a
!> spherical basis R(r)·Y_l^m(θ,φ). The grid stores composite indices
!> (r, l, m), and the InitFunction computes the radial part for matching
!> (l, m) quantum numbers.
!>
!> Available Sub-Backends
!> ----------------------
!> - **HydrogenLike** : Hydrogenic radial functions R_nl(r) for Coulomb potentials.
!>                      Quantum numbers (n, l, m) specified per orbital via JSON.
!>
!> Procedure Pointer Contract
!> --------------------------
!> `OrbsInit_Ylm_InitFunction(r, l, m, index, bt_)` returns complex amplitude
!> at radius r for the given angular quantum numbers (l, m) and orbital index.
!> Returns zero if (l, m) doesn't match the orbital's target quantum numbers.
!>
!> JSON Configuration
!> ------------------
!>   "orbsInit": {
!>     "ylm": {
!>       "hydrogenLike": {
!>         "charge": 1.0,           // effective nuclear charge Z
!>         "n": [1, 2, 2],          // principal quantum numbers
!>         "l": [0, 0, 1],          // orbital angular momentum
!>         "m": [0, 0, 0]           // magnetic quantum numbers
!>       }
!>     }
!>   }
!>
!> @note Grid coordinates accessed via Grid_Ylm_rCoord, Grid_Ylm_lCoord, Grid_Ylm_mCoord.
!> @see M_OrbsInit_Ylm_HydrogenLike for the hydrogenic radial implementation.
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
