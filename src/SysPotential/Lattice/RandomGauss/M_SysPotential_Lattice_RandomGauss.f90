! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> RandomGauss potential on lattice grids.
!>
!> Behaviour:
!> - Time-independent: V(i) = fixed randomGauss values drawn once from the chosen
!>   distribution during setup and cached for all subsequent calls
module M_SysPotential_Lattice_RandomGauss
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Parse parameters and wire the lattice randomGauss implementation.
    !>
    !> Reads the configuration (typically from JSON), initializes the module
    !> data (distribution parameters), and connects the implementation into
    !> the lattice backend and grid-agnostic pointers.
    module subroutine SysPotential_Lattice_RandomGauss_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Gaussian distribution mean value.
  real(R64) :: SysPotential_Lattice_RandomGauss_meanValue = 0.0_R64
  !> Gaussian distribution standard deviation.
  real(R64) :: SysPotential_Lattice_RandomGauss_stdValue = 1.0_R64
  !> Seed for the random number generator. If -1, the seed is generated from the system clock.
  integer(I32) :: SysPotential_Lattice_RandomGauss_seed = -1
  !> Cached random site values generated in Setup
  real(R64), allocatable :: SysPotential_Lattice_RandomGauss_randomGaussValues(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
