! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Linear-grid absorber family and factory.
!>
!> This module collects absorber implementations that operate on a 1D linear
!> grid (e.g., cosine masks) and exposes a factory routine that selects and
!> installs one of them into the public interface `M_Absorber`.
module M_Absorber_Linear
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Install a linear-grid absorber implementation and read input.
    !>
    !> Binds the procedure pointers in `M_Absorber` to a specific linear-grid
    !> implementation (for example, a cosine profile) and consumes user
    !> configuration (e.g., from JSON). Call this during program initialization.
    module subroutine Absorber_Linear_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
