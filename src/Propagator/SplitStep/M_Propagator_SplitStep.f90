! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Split-operator propagator backend (framework for higher-order schemes).
!>
!> This backend decomposes the Hamiltonian into parts (e.g., \(\hat{H}=\hat{A}+\hat{B}\))
!> that are each cheap to apply (e.g., diagonal in complementary bases) and
!> composes short-time evolutions using Suzuki–Trotter factorizations.
!>
!> Example (second order):
!> \[
!>   e^{-i\hat{H}\Delta t} \approx e^{-i\hat{A}\Delta t/2}\, e^{-i\hat{B}\Delta t}\, e^{-i\hat{A}\Delta t/2} + \mathcal{O}(\Delta t^3)
!> \]
!>
!> Specific composition orders are implemented in submodules such as
!> `M_Propagator_SplitStep_Order2` and `M_Propagator_SplitStep_Order4`.
module M_Propagator_SplitStep
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate and wire the split-operator backend.
    !>
    !> Assigns `M_Propagator` procedure pointers to split-step implementations
    !> and configures the selected composition order (via this module or
    !> order-specific submodules). After fabrication, call
    !> `M_Propagator::Propagator_Setup` and then `Propagator_Propagate`.
    module subroutine Propagator_SplitStep_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
