! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Second-order split-operator propagator backend.
!>
!> Provides a symmetric Suzuki–Trotter composition with local error
!> \(\mathcal{O}(\Delta t^3)\) and global error \(\mathcal{O}(\Delta t^2)\):
!> \[
!>   e^{-i\hat{H}\Delta t} \approx e^{-i\hat{A}\Delta t/2}\, e^{-i\hat{B}\Delta t}\, e^{-i\hat{A}\Delta t/2}.
!> \]
!> This order balances simplicity and accuracy and is often a robust default
!> when \(\hat{A}\) and \(\hat{B}\) are efficiently applicable.
module M_Propagator_SplitStep_Order2
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate and wire the order-2 split-operator backend.
    !>
    !> Assigns `M_Propagator` procedure pointers to the second-order
    !> composition implementation. After fabrication, call setup and then
    !> propagate via the `M_Propagator` facade.
    module subroutine Propagator_SplitStep_Order2_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
