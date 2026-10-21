! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Adams–Lawson predictor–corrector scheme of the exponential integrators
!> (exponential Lawson multistep method, Auzinger et al., J. Comput. Math.
!> Data Sci. 1, 100014 (2021)).
!>
!> L is removed by the integrating factor E(τ) = exp(τ L): in the
!> interaction picture v = E(tRef - t) u the remaining derivative is
!>   v' = G(t, v) = E(tRef - t) N(t, E(t - tRef) v),
!> integrated by a k-step Adams–Bashforth predictor and Adams–Moulton
!> corrector (PECE, order k + 1). Per step: two N and four exact
!> propagations; the step size must stay fixed (a change restarts the
!> history).
!>
!> G carries the phases of L, so the polynomial fit of the history limits the
!> step size like for the exponential Adams scheme; on the Be HHG benchmark
!> it stayed stable at a spectral cut where `adams` was not, but was less
!> accurate.
module M_Propagator_Exponential_AdamsLawson
  use M_Utils_Types

  implicit none

  interface
    !> Reads `"propagator.exponential.adamsLawson.nSteps"` and binds
    !> `Propagator_Propagate`.
    module subroutine Propagator_Exponential_AdamsLawson_Fabricate
    end subroutine
  end interface

end module
