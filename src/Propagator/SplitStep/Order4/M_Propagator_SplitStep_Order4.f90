! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fourth-order split-operator propagator backend.
!>
!> Provides a higher-order symmetric composition (e.g., Yoshida/Forest–Ruth
!> type) with global error \(\mathcal{O}(\Delta t^4)\). Coefficients are
!> chosen to cancel lower-order error terms while maintaining time symmetry.
!> A representative structure is a weighted sequence of half-steps of \(\hat{A}\)
!> and full steps of \(\hat{B}\):
!> \[
!>   U_4(\Delta t) = \prod_k e^{-i a_k \hat{A}\, \Delta t}\, e^{-i b_k \hat{B}\, \Delta t}.
!> \]
!> This backend is well-suited when higher accuracy per step offsets the cost
!> of additional operator applications.
module M_Propagator_SplitStep_Order4
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate and wire the order-4 split-operator backend.
    !>
    !> Assigns `M_Propagator` procedure pointers to the fourth-order composition
    !> implementation. After fabrication, call setup and then propagate via the
    !> `M_Propagator` facade.
    module subroutine Propagator_SplitStep_Order4_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
