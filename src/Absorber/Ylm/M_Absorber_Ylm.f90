! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Ylm-grid absorber family and factory.
!>
!> @details
!> This module collects absorber implementations that operate on a Ylm grid
!> (radial grid with spherical-harmonic (l,m) channels) and exposes a factory
!> routine that selects and installs one of them into the public interface
!> `M_Absorber`.
!>
!> Ylm absorbers apply a purely radial mask M(r): since the flattened Ylm grid
!> stores the radius of every (r,l,m) point, the same radial mask acts on each
!> (l,m) channel independently.
!>
!> @par Available Implementations
!> - **cosinus**: Smooth cosine^(1/n) radial taper (see `M_Absorber_Ylm_Cosinus`)
!>
!> @see M_Absorber, M_Absorber_Ylm_Cosinus
module M_Absorber_Ylm
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Install a Ylm-grid absorber implementation and read input.
    !>
    !> @details
    !> Binds the procedure pointers in `M_Absorber` to a specific Ylm-grid
    !> implementation based on JSON configuration. Currently dispatches to:
    !>
    !> - **absorber.ylm.cosinus**: Cosine-profile radial absorber
    !>
    !> @pre  JSON configuration must specify one of the supported variants.
    !> @post `Absorber_Setup` and `Absorber_ApplyAbsorber` are bound.
    module subroutine Absorber_Ylm_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
