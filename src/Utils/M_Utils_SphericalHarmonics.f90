! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Spherical harmonics helpers and Gaunt coefficient evaluation.
!>
!> Provides utilities to compute integrals involving products of spherical
!> harmonics via Wigner 3j symbols and related special functions.
module M_Utils_SphericalHarmonics
  use M_Utils_Types  ! Added import for R64 and I32 types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pure function SphericalHarmonics_GauntCoefficient(l1, m1, l2, m2, l3, m3) result(res)
    use M_Utils_Constants
    use M_Utils_Constants
    use M_Utils_SfGslLib
    ! Computes the Gaunt coefficient, which is the integral of three spherical harmonics:
    !
    ! ∫₀²ᵖ ∫₀ᵖ Yₗ₁,ₘ₁(θ,ϕ) Yₗ₂,ₘ₂(θ,ϕ) Y^*ₗ₃,ₘ₃(θ,ϕ) sin(θ) dθ dϕ
    !
    ! The formula uses Wigner 3j symbols:
    ! √[(2l₁+1)(2l₂+1)(2l₃+1)/(4π)] * (l₁ l₂ l₃; 0 0 0) * (l₁ l₂ l₃; m₁ m₂ m₃)

    real(R64) :: res
    integer(I32), intent(in) :: l1, m1, l2, m2, l3, m3

    real(R64) :: prefactor, iTmp, jTmp
    integer(I32) :: phase

    ! Apply phase for complex conjugation
    if (mod(abs(m3), 2) .eq. 0) then
      phase = 1
    else
      phase = -1
    end if

    ! Calculate the prefactor
    prefactor = phase * sqrt((2.0_R64 * l1 + 1.0_R64) * (2.0_R64 * l2 + 1.0_R64) * &
                             (2.0_R64 * l3 + 1.0_R64) / (4.0_R64 * PI))

    ! Calculate the first Wigner 3j symbol (l₁ l₂ l₃; 0 0 0)
    iTmp = SfGslLib_Wigner3jGsl(l1, l2, l3, 0, 0, 0)

    ! Calculate the second Wigner 3j symbol (l₁ l₂ l₃; m₁ m₂ m₃)
    jTmp = SfGslLib_Wigner3jGsl(l1, l2, l3, m1, m2, -m3) ! Note the sign because y3 is conjugated

    ! Calculate the Gaunt coefficient
    res = prefactor * iTmp * jTmp

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module
