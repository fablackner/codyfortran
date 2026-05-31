! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the linear-grid absorber fabrication logic.
submodule(M_Absorber_Linear) S_Absorber_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Factory routine for linear-grid absorber variants.
  !>
  !> @details
  !> Inspects JSON configuration and delegates to the appropriate variant:
  !>
  !> - **absorber.linear.cosinus**: Cosine-profile mask
  !>
  !> Terminates with error if no valid variant is specified.
  module subroutine Absorber_Linear_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Absorber_Linear_Cosinus

    call Say_Fabricate("absorber.linear")

    !------------------------------------
    ! branch: select linear absorber variant
    !------------------------------------

    if (Json_GetExistence("absorber.linear.cosinus")) then
      call Absorber_Linear_Cosinus_Fabricate

    else
      error stop "absorber.linear is missing one of: cosinus"
    end if

  end subroutine

end submodule
