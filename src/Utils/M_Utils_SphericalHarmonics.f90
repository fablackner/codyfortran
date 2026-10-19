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
  !> (GauntIndex(l1,m1), GauntIndex(l2,m2), l3+1) with
  !> SphericalHarmonics_GauntIndex over |m1|, |m2| <= the tabulated m bound. Read-only after setup, so lookups are thread-safe;
  !> building it concurrently with lookups is not.
  real(R64), allocatable :: SphericalHarmonics_gauntTable(:, :, :)

  !> Tabulated bounds; -1 while no table exists
  integer(I32) :: SphericalHarmonics_gauntTableL1 = -1
  integer(I32) :: SphericalHarmonics_gauntTableL2 = -1
  integer(I32) :: SphericalHarmonics_gauntTableL3 = -1
  integer(I32) :: SphericalHarmonics_gauntTableM = -1

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Builds (or grows) the Gaunt coefficient table so that all combinations
  !> with l1 <= l1max, l2 <= l2max, l3 <= l3max and |m1|, |m2| <= mmax_
  !> (default: no m restriction) are tabulated. Call during
  !> fabrication/setup (serial context) before any tabulated lookups; repeated
  !> calls only rebuild when a bound grows.
  subroutine SphericalHarmonics_EnsureGauntTable(l1max, l2max, l3max, mmax_)

    integer(I32), intent(in) :: l1max
    integer(I32), intent(in) :: l2max
    integer(I32), intent(in) :: l3max
    integer(I32), intent(in), optional :: mmax_

    integer(I32) :: newL1, newL2, newL3, newM, mmax
    integer(I32) :: l1, m1, l2, m2, l3
    integer(I32) :: i1, i2

    mmax = max(l1max, l2max)
    if (present(mmax_)) mmax = min(mmax_, mmax)

    newL1 = max(l1max, SphericalHarmonics_gauntTableL1)
    newL2 = max(l2max, SphericalHarmonics_gauntTableL2)
    newL3 = max(l3max, SphericalHarmonics_gauntTableL3)
    newM = max(mmax, SphericalHarmonics_gauntTableM)

    ! Already covered
    if (allocated(SphericalHarmonics_gauntTable) .and. &
        newL1 <= SphericalHarmonics_gauntTableL1 .and. &
        newL2 <= SphericalHarmonics_gauntTableL2 .and. &
        newL3 <= SphericalHarmonics_gauntTableL3 .and. &
        newM <= SphericalHarmonics_gauntTableM) return

    if (allocated(SphericalHarmonics_gauntTable)) deallocate (SphericalHarmonics_gauntTable)
    allocate (SphericalHarmonics_gauntTable(SphericalHarmonics_GauntIndex(newL1 + 1, -min(newL1 + 1, newM), newM) - 1, &
                                            SphericalHarmonics_GauntIndex(newL2 + 1, -min(newL2 + 1, newM), newM) - 1, &
                                            newL3 + 1))
    SphericalHarmonics_gauntTable = 0.0_R64

    !$omp parallel do default(shared) private(l1, m1, l2, m2, l3, i1, i2) schedule(dynamic)
    do l1 = 0, newL1
      do m1 = -min(l1, newM), min(l1, newM)
        i1 = SphericalHarmonics_GauntIndex(l1, m1, newM)
        do l2 = 0, newL2
          do m2 = -min(l2, newM), min(l2, newM)
            i2 = SphericalHarmonics_GauntIndex(l2, m2, newM)
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
    SphericalHarmonics_gauntTableM = newM

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Compact table index of (l, m) for |m| <= min(l, mmax): channels are
  !> ordered by l, then m, and each l contributes 2 min(l, mmax) + 1 entries
  !> (l² + l + m + 1 when mmax >= l).
  pure function SphericalHarmonics_GauntIndex(l, m, mmax) result(res)

    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    integer(I32), intent(in) :: mmax
    integer(I32) :: res

    if (l <= mmax + 1) then
      res = l * l
    else
      res = (mmax + 1)**2 + (l - mmax - 1) * (2 * mmax + 1)
    end if
    res = res + m + min(l, mmax) + 1

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Gaunt coefficient for the (only nonvanishing) case m3 = m1 + m2, served
  !> from the precomputed table when in range, otherwise computed directly.
  pure function SphericalHarmonics_GauntTabulated(l1, m1, l2, m2, l3) result(res)

    real(R64) :: res
    integer(I32), intent(in) :: l1, m1, l2, m2, l3

    if (l1 <= SphericalHarmonics_gauntTableL1 .and. &
        l2 <= SphericalHarmonics_gauntTableL2 .and. &
        l3 <= SphericalHarmonics_gauntTableL3 .and. &
        abs(m1) <= SphericalHarmonics_gauntTableM .and. &
        abs(m2) <= SphericalHarmonics_gauntTableM) then
      res = SphericalHarmonics_gauntTable(SphericalHarmonics_GauntIndex(l1, m1, SphericalHarmonics_gauntTableM), &
                                          SphericalHarmonics_GauntIndex(l2, m2, SphericalHarmonics_gauntTableM), l3 + 1)
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
