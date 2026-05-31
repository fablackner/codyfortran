! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for hydrogen-like radial orbital initialization.
submodule(M_OrbsInit_Ylm_HydrogenLike) S_OrbsInit_Ylm_HydrogenLike

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Read quantum number arrays from JSON and bind the InitFunction pointer.
  module subroutine OrbsInit_Ylm_HydrogenLike_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit_Ylm

    implicit none

    real(R64) :: charge

    call Say_Fabricate("orbsInit.ylm.hydrogenLike")

    ! Read effective nuclear charge (default: hydrogen Z=1)
    charge = Json_Get("orbsInit.ylm.hydrogenLike.charge", 1.0_R64)
    OrbsInit_Ylm_HydrogenLike_charge = charge

    ! Read quantum number arrays (one entry per orbital)
    OrbsInit_Ylm_HydrogenLike_n = Json_Get("orbsInit.ylm.hydrogenLike.n", [1])
    OrbsInit_Ylm_HydrogenLike_l = Json_Get("orbsInit.ylm.hydrogenLike.l", [0])
    OrbsInit_Ylm_HydrogenLike_m = Json_Get("orbsInit.ylm.hydrogenLike.m", [0])

    OrbsInit_Ylm_InitFunction => InitFunction

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Compute hydrogen-like radial wavefunction at radius r.
!>
!> @details
!> Evaluates R_nl(r) = ρ^l · L_{n-l-1}^{2l+1}(ρ) · exp(-ρ/2) where ρ = 2Zr/(na₀).
!> Returns zero if the grid point's (l, m) doesn't match the orbital's quantum numbers.
!>
!> @param[in]  r      Radial coordinate
!> @param[in]  l      Angular momentum from grid point
!> @param[in]  m      Magnetic quantum number from grid point
!> @param[in]  index  Orbital index (selects quantum numbers from config arrays)
!> @param[in]  bt_    Body type (unused, for interface compatibility)
!> @return     res    Radial wavefunction value (zero if l,m don't match)
  function InitFunction(r, l, m, index, bt_) result(res)
    use M_Utils_UnusedVariables
    use M_Utils_SfGslLib

    complex(R64)             :: res
    real(R64), intent(in)    :: r
    integer(I32), intent(in) :: l, m, index
    integer(I32), intent(in), optional :: bt_

    real(R64)    :: charge, a0, rho
    integer(I32) :: n

    if (.false.) call UnusedVariables_Mark(m, bt_)

    res = (0.0_R64, 0.0_R64)

    ! Check if (l, m) matches this orbital's quantum numbers
    if (l /= OrbsInit_Ylm_HydrogenLike_l(index)) return
    if (m /= OrbsInit_Ylm_HydrogenLike_m(index)) return

    n = OrbsInit_Ylm_HydrogenLike_n(index)
    charge = OrbsInit_Ylm_HydrogenLike_charge
    a0 = 1.0_R64  ! Bohr radius in atomic units

    ! Dimensionless radial coordinate ρ = 2Zr/(na₀)
    rho = 2.0_R64 * charge * r / (n * a0)

    ! Radial wavefunction: R_nl(r) ∝ ρ^l · L_{n-l-1}^{2l+1}(ρ) · exp(-ρ/2)
    res = exp(-rho / 2.0_R64) * (rho**l) * &
          SfGslLib_Laguerre(n - l - 1, 2 * l + 1, rho)

  end function

end submodule
