! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the FEM-grid absorber fabrication logic.
submodule(M_Absorber_Fem) S_Absorber_Fem

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Factory routine for FEM-grid absorber variants.
  !>
  !> @details
  !> - **absorber.fem.cosinus**: cosine-profile boundary layers
  !>
  !> Terminates with an error if grid.fem is missing or no valid variant is
  !> specified.
  module subroutine Absorber_Fem_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Absorber_Fem_Cosinus

    call Say_Fabricate("absorber.fem")

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.fem")) then
      error stop "grid.fem is required for absorber.fem"
    end if

    !------------------------------------
    ! branch: select fem absorber variant
    !------------------------------------

    if (Json_GetExistence("absorber.fem.cosinus")) then
      call Absorber_Fem_Cosinus_Fabricate

    else
      error stop "absorber.fem is missing one of: cosinus"
    end if

  end subroutine

end submodule
