! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the absorber fabrication logic.
!>
!> @details
!> This submodule reads the JSON configuration and dispatches to the
!> appropriate absorber family (e.g., Linear). If no absorber section is
!> present, a no-op absorber is installed.
submodule(M_Absorber) S_Absorber

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Factory routine that selects and installs an absorber implementation.
  !>
  !> @details
  !> Inspects the JSON configuration tree and delegates to the appropriate
  !> sub-factory. Currently supported branches:
  !>
  !> - **absorber.linear**: 1D linear-grid absorbers (e.g., cosine profile)
  !> - **absorber.ylm**: Ylm-grid radial absorbers (e.g., cosine profile)
  !>
  !> If no absorber block is found, the no-op fallback is installed, leaving
  !> wavefunctions unmodified. This is useful for testing or when boundary
  !> reflections are not a concern.
  module subroutine Absorber_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Absorber_Linear
    use M_Absorber_Ylm

    call Say_Fabricate("absorber")

    !------------------------------------
    ! branch: select absorber family
    !------------------------------------

    if (Json_GetExistence("absorber.linear")) then
      call Absorber_Linear_Fabricate

    else if (Json_GetExistence("absorber.ylm")) then
      call Absorber_Ylm_Fabricate

    else
      ! No absorber configured → install no-op (wavefunctions pass through)
      write (*, '(1X, A, A, A)') red, 'de: Absorber => noAbsorber', reset
      Absorber_ApplyAbsorber => NoOpProcedures_ApplyAbsorber

    end if

  end subroutine

end submodule
