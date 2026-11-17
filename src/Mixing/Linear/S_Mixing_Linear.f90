! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Linear mixing implementation submodule.
submodule(M_Mixing_Linear) S_Mixing_Linear

  implicit none

  !> Damping parameter (0 < α ≤ 1). Read from JSON during fabrication.
  real(R64) :: alpha = 1.0_R64

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Mixing_Linear_Fabricate()
    use M_Utils_Say
    use M_Utils_Json
    use M_Mixing

    call Say_Fabricate("mixing.linear")

    !------------------------------------
    ! read configuration parameters
    !------------------------------------

    alpha = Json_Get("alpha", 1.0_R64, path_="mixing.linear")
    if (.not. (alpha > 0.0_R64 .and. alpha <= 1.0_R64)) error stop "mixing.linear.alpha must be in (0, 1]"

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Mixing_Mix => Mix
    Mixing_Reset => ResetMixer

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Damped linear step: xMixed = (1-α)·xMixed + α·xRaw
  subroutine Mix(xMixed, xRaw)
    complex(R64), intent(inout), contiguous :: xMixed(:)
    complex(R64), intent(in), contiguous :: xRaw(:)

    xMixed = (1.0_R64 - alpha) * xMixed + alpha * xRaw

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Stateless: nothing to reset.
  subroutine ResetMixer
  end subroutine

end submodule
