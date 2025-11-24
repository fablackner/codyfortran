! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module M_OrbsInit_Ylm_HydrogenLike
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate a hydrogen-like (Coulombic) Ylm initializer.
    !> Configures radial part parameters for hydrogen-like orbitals and wires the
    !> callouts in the Ylm backend. Quantum numbers may be provided per-orbital.
    module subroutine OrbsInit_Ylm_HydrogenLike_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Effective nuclear charge Z used in the hydrogen-like radial model.
  real(R64) :: OrbsInit_Ylm_HydrogenLike_charge
  !> Principal quantum number n for each orbital to initialize.
  integer(I32), allocatable :: OrbsInit_Ylm_HydrogenLike_n(:)
  !> Orbital angular momentum quantum number l for each orbital.
  integer(I32), allocatable :: OrbsInit_Ylm_HydrogenLike_l(:)
  !> Magnetic quantum number m for each orbital.
  integer(I32), allocatable :: OrbsInit_Ylm_HydrogenLike_m(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
