! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the Ylm-grid absorber fabrication logic.
submodule(M_Absorber_Ylm) S_Absorber_Ylm

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Factory routine for Ylm-grid absorber variants.
  !>
  !> @details
  !> Inspects JSON configuration and delegates to the appropriate variant:
  !>
  !> - **absorber.ylm.cosinus**: Cosine-profile radial mask
  !>
  !> Terminates with error if no valid variant is specified.
  module subroutine Absorber_Ylm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Absorber_Ylm_Cosinus

    call Say_Fabricate("absorber.ylm")

    !------------------------------------
    ! branch: select ylm absorber variant
    !------------------------------------

    if (Json_GetExistence("absorber.ylm.cosinus")) then
      call Absorber_Ylm_Cosinus_Fabricate

    else
      error stop "absorber.ylm is missing one of: cosinus"
    end if

  end subroutine

end submodule
