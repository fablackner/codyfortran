! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Generic (spin-unresolved) Valdemoro reconstruction scheme.
!>
!> @details
!> Implements the Valdemoro approximation for the full 2-RDM without spin symmetry.
!> The 3-RDM is reconstructed by setting the connected 3-body cumulant to zero (C³=0),
!> yielding the "free-particle" approximation:
!>
!>     D³(123) ≈ A[D¹(1)D¹(2)D¹(3)] + A[D¹(1)Δ²(23)] + A[Δ²(12)D¹(3)]
!>
!> where A denotes antisymmetrization and Δ² = D² - A[D¹D¹] is the 2-body cumulant.
!>
!> **JSON path**: `reconstruction.generic.genericValdemoro`
!>
!> @see C. Valdemoro, Phys. Rev. A 45, 4462 (1992)
module M_Reconstruction_Generic_Valdemoro
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the generic Valdemoro implementation at runtime.
    !>
    !> Wires `Reconstruction_Generic_FillCollisionOp2` in
    !> `M_Reconstruction_Generic` to the Valdemoro-specific routine and may
    !> perform light initialization specific to this variant.
    module subroutine Reconstruction_Generic_Valdemoro_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! local data
  !=============================================================================

end module
