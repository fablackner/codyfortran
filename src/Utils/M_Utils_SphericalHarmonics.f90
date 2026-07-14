! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Spherical harmonics helpers and Gaunt coefficient evaluation.
!>
!> Provides utilities to compute integrals involving products of spherical
!> harmonics via Wigner 3j symbols and related special functions.
!>
!> A precomputed Gaunt coefficient table (built once at fabrication via
!> SphericalHarmonics_EnsureGauntTable) eliminates the per-call GSL Wigner-3j
!> evaluations from hot loops: SphericalHarmonics_GauntTabulated reads the
!> table and falls back to the direct computation outside the tabulated range.
module M_Utils_SphericalHarmonics
  use M_Utils_Types  ! Added import for R64 and I32 types

  implicit none

  !> Gaunt table for m3 = m1 + m2 (the only nonvanishing case), indexed as
  !> (l1²+l1+m1+1, l2²+l2+m2+1, l3+1). Read-only after setup, so lookups are
  !> thread-safe; building it concurrently with lookups is not.
  real(R64), allocatable :: SphericalHarmonics_gauntTable(:, :, :)

  !> Tabulated bounds; -1 while no table exists
  integer(I32) :: SphericalHarmonics_gauntTableL1 = -1
  integer(I32) :: SphericalHarmonics_gauntTableL2 = -1
  integer(I32) :: SphericalHarmonics_gauntTableL3 = -1

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Builds (or grows) the Gaunt coefficient table so that all combinations
  !> with l1 <= l1max, l2 <= l2max, l3 <= l3max are tabulated. Call during
  !> fabrication/setup (serial context) before any tabulated lookups; repeated
  !> calls only rebuild when a bound grows.
  subroutine SphericalHarmonics_EnsureGauntTable(l1max, l2max, l3max)

    integer(I32), intent(in) :: l1max
    integer(I32), intent(in) :: l2max
    integer(I32), intent(in) :: l3max

    integer(I32) :: newL1, newL2, newL3
    integer(I32) :: l1, m1, l2, m2, l3
    integer(I32) :: i1, i2

    newL1 = max(l1max, SphericalHarmonics_gauntTableL1)
    newL2 = max(l2max, SphericalHarmonics_gauntTableL2)
    newL3 = max(l3max, SphericalHarmonics_gauntTableL3)

    ! Already covered
    if (allocated(SphericalHarmonics_gauntTable) .and. &
        newL1 <= SphericalHarmonics_gauntTableL1 .and. &
        newL2 <= SphericalHarmonics_gauntTableL2 .and. &
        newL3 <= SphericalHarmonics_gauntTableL3) return

    if (allocated(SphericalHarmonics_gauntTable)) deallocate (SphericalHarmonics_gauntTable)
    allocate (SphericalHarmonics_gauntTable((newL1 + 1)**2, (newL2 + 1)**2, newL3 + 1))
    SphericalHarmonics_gauntTable = 0.0_R64

    !$omp parallel do default(shared) private(l1, m1, l2, m2, l3, i1, i2) schedule(dynamic)
    do l1 = 0, newL1
      do m1 = -l1, l1
        i1 = l1 * l1 + l1 + m1 + 1
        do l2 = 0, newL2
          do m2 = -l2, l2
            i2 = l2 * l2 + l2 + m2 + 1
            do l3 = abs(l1 - l2), min(l1 + l2, newL3)
              ! Selection rules: parity and |m3| <= l3 (with m3 = m1 + m2)
              if (mod(l1 + l2 + l3, 2) .ne. 0) cycle
              if (abs(m1 + m2) > l3) cycle
              SphericalHarmonics_gauntTable(i1, i2, l3 + 1) = &
                SphericalHarmonics_GauntCoefficient(l1, m1, l2, m2, l3, m1 + m2)
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do

    SphericalHarmonics_gauntTableL1 = newL1
    SphericalHarmonics_gauntTableL2 = newL2
    SphericalHarmonics_gauntTableL3 = newL3

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Gaunt coefficient for the (only nonvanishing) case m3 = m1 + m2, served
  !> from the precomputed table when in range, otherwise computed directly.
  pure function SphericalHarmonics_GauntTabulated(l1, m1, l2, m2, l3) result(res)

    real(R64) :: res
    integer(I32), intent(in) :: l1, m1, l2, m2, l3

    if (l1 <= SphericalHarmonics_gauntTableL1 .and. &
        l2 <= SphericalHarmonics_gauntTableL2 .and. &
        l3 <= SphericalHarmonics_gauntTableL3) then
      res = SphericalHarmonics_gauntTable(l1 * l1 + l1 + m1 + 1, l2 * l2 + l2 + m2 + 1, l3 + 1)
    else
      res = SphericalHarmonics_GauntCoefficient(l1, m1, l2, m2, l3, m1 + m2)
    end if

  end function

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
