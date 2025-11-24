! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Full Configuration Interaction (Full CI) with static orbitals.
!>
!> Represents the many-body state with the complete set of Slater determinants
!> and time-dependent CI coefficients, while keeping one-body orbitals fixed.
!> Provides the exact TDSE solution within the chosen orbital basis.
!>
!> State layout
!> - The packed state contains CI coefficients ordered consistently with the
!>   determinant/configuration list implied by the selected basis and body types.
module M_Method_Mb_OrbBased_Tdci
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds Full CI procedure pointers at runtime and initializes the CI state.
    !> Typically sets `Method_Setup`, `Method_GetEnergy`, and `Method_TimeDerivative`.
    module subroutine Method_Mb_OrbBased_Tdci_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
