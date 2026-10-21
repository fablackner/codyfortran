! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Exponential Adams predictor–corrector scheme of the exponential
!> integrators (Hochbruck & Ostermann, Acta Numerica 19, 209 (2010), Sec. 2.5).
!>
!> The variation-of-constants formula
!>   u(t_n + h) = exp(h L) u_n + ∫_0^h exp((h - s) L) N(t_n + s) ds
!> is integrated with N replaced by its interpolation polynomial through the
!> last k values (predictor) and through the predicted value as well
!> (corrector). Per step: two N and two phi combinations; the step size must
!> stay fixed (a change restarts the history).
!>
!> N is interpolated in time, so components of N that oscillate fast with the
!> phases of L (e.g. high-energy orbital components coupled by N) limit the
!> step size: h times the largest exactly propagated frequency should stay
!> well below 1 (for the spectral cut: h * method.mb.orbBased.exactFieldFreeCut).
module M_Propagator_Exponential_Adams
  use M_Utils_Types

  implicit none

  interface
    !> Reads `"propagator.exponential.adams.nSteps"` and binds
    !> `Propagator_Propagate`.
    module subroutine Propagator_Exponential_Adams_Fabricate
    end subroutine
  end interface

end module
