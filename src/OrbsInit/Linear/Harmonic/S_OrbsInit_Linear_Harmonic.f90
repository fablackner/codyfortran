! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_OrbsInit_Linear_Harmonic) S_OrbsInit_Linear_Harmonic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine OrbsInit_Linear_Harmonic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit
    use M_OrbsInit_Linear

    implicit none

    real(R64) :: position
    real(R64) :: omega

    call Say_Fabricate("orbsInit.linear.harmonic")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    position = Json_Get("orbsInit.linear.harmonic.position", 0.0_R64)
    omega = Json_Get("orbsInit.linear.harmonic.omega", 1.0_R64)
    OrbsInit_Linear_Harmonic_position = position
    OrbsInit_Linear_Harmonic_omega = omega

    OrbsInit_Linear_InitFunction => InitFunction

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function InitFunction(x, index, bt_) result(res)
    use M_Utils_UnusedVariables
    use M_Utils_SfGslLib

    real(R64)                :: res
    real(R64), intent(in)    :: x
    integer(I32), intent(in) :: index
    integer(I32), intent(in), optional :: bt_

    real(R64) :: position, omega

    integer(I32) :: nx
    real(R64)    :: a

    if (.false.) call UnusedVariables_Mark(bt_)

    ! Set quantum number directly based on index
    nx = index - 1

    ! Set parameters for harmonic oscillator
    position = OrbsInit_Linear_Harmonic_position
    omega = OrbsInit_Linear_Harmonic_omega
    a = 1.0_R64 / sqrt(omega)

    res = exp(-0.5_R64 * ((x - position) / a)**2) * SfGslLib_Hermite((x - position) / a, nx)

  end function

end submodule
