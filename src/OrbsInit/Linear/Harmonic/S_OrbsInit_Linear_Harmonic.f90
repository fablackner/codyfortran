! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for harmonic oscillator orbital initialization.
submodule(M_OrbsInit_Linear_Harmonic) S_OrbsInit_Linear_Harmonic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Read harmonic parameters from JSON and bind the InitFunction pointer.
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
!> @brief Compute harmonic oscillator eigenstate amplitude at position x.
!>
!> @details
!> Evaluates the (unnormalized) harmonic oscillator wavefunction:
!>   ψ_n(x) = H_n(ξ) · exp(-ξ²/2)
!> where ξ = (x - x₀)/a, a = 1/√ω, and n = index - 1.
!>
!> @param[in]  x      Spatial coordinate
!> @param[in]  index  Orbital index (1-based; n = index - 1)
!> @param[in]  bt_    Body type (unused, for interface compatibility)
!> @return     res    Unnormalized wavefunction amplitude (real-valued)
  function InitFunction(x, index, bt_) result(res)
    use M_Utils_UnusedVariables
    use M_Utils_SfGslLib

    real(R64)                :: res
    real(R64), intent(in)    :: x
    integer(I32), intent(in) :: index
    integer(I32), intent(in), optional :: bt_

    real(R64)    :: position, omega, a, xi
    integer(I32) :: n

    if (.false.) call UnusedVariables_Mark(bt_)

    ! Map 1-based index to quantum number n (ground state n=0 at index=1)
    n = index - 1

    ! Retrieve oscillator parameters
    position = OrbsInit_Linear_Harmonic_position
    omega = OrbsInit_Linear_Harmonic_omega
    a = 1.0_R64 / sqrt(omega)  ! oscillator length scale

    ! Dimensionless coordinate
    xi = (x - position) / a

    ! Harmonic oscillator eigenstate: H_n(ξ) · exp(-ξ²/2)
    res = exp(-0.5_R64 * xi**2) * SfGslLib_Hermite(xi, n)

  end function

end submodule
