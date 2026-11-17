! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Krogstad scheme of the exponential integrators: exponential Runge–Kutta of
!> order 4 (S. Krogstad, J. Comput. Phys. 203, 72 (2005)), one step with four
!> evaluations of N.
!>
!> The step itself is the building block `Propagator_Exponential_KrogstadStep`
!> of M_Propagator_Exponential, shared with the multistep schemes, which start
!> their history with it.
module M_Propagator_Exponential_Krogstad4
  use M_Utils_Types

  implicit none

  interface
    !> Binds `Propagator_Propagate` to the Krogstad step
    !> (`"propagator.exponential.krogstad4"` has no keys).
    module subroutine Propagator_Exponential_Krogstad4_Fabricate
    end subroutine
  end interface

end module
