! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Cosine-profile absorber for a linear grid.
!>
!> This module implements a smooth, reflection-minimizing mask based on a
!> cosine profile for both ends of a 1D linear grid. The mask is unity in the
!> interior region and decays smoothly to zero inside boundary layers of a
!> configurable thickness, thereby absorbing wavefunctions and reducing
!> spurious reflections.
!>
!> Conceptually, the mask can be expressed as
!> $M(x) = \cos^n\!\left(\tfrac{\pi}{2}\,\xi(x)\right)$ for $0\le \xi \le 1$,
!> with $\xi(x)$ a linear mapping from the start of the absorbing layer to the
!> boundary, and $n$ an integer order controlling the taper steepness. Outside
!> the absorbing layers $M(x)=1$ and at the boundaries $M(x)=0$.
!>
!> This module only installs the implementation by wiring procedure pointers in
!> the parent `M_Absorber` module; it does not expose a new public API.
module M_Absorber_Linear_Cosinus
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Install the cosine-profile linear absorber and read its configuration.
    !>
    !> Binds `M_Absorber`'s procedure pointers to this cosine implementation and
    !> reads user-provided parameters such as the absorbing-layer thickness and
    !> the power/order of the cosine taper (n).
    module subroutine Absorber_Linear_Cosinus_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
