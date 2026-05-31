! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Linear-grid absorber family and factory.
!>
!> @details
!> This module collects absorber implementations that operate on a 1D linear
!> grid (e.g., cosine masks) and exposes a factory routine that selects and
!> installs one of them into the public interface `M_Absorber`.
!>
!> Linear absorbers assume a symmetric domain [-L, +L] and apply damping
!> in regions where |x| exceeds a configurable onset threshold.
!>
!> @par Available Implementations
!> - **cosinus**: Smooth cosine^(1/n) taper (see `M_Absorber_Linear_Cosinus`)
!>
!> @see M_Absorber, M_Absorber_Linear_Cosinus
module M_Absorber_Linear
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Install a linear-grid absorber implementation and read input.
    !>
    !> @details
    !> Binds the procedure pointers in `M_Absorber` to a specific linear-grid
    !> implementation based on JSON configuration. Currently dispatches to:
    !>
    !> - **absorber.linear.cosinus**: Cosine-profile absorber
    !>
    !> @pre  JSON configuration must specify one of the supported variants.
    !> @post `Absorber_Setup` and `Absorber_ApplyAbsorber` are bound.
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
