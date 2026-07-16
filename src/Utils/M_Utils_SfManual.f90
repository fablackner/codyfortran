! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Manual implementations of selected special functions.
module M_Utils_SfManual
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure recursive function SfManual_Factorial(n) result(res)

    integer(I32) :: res
    integer(I32), intent(in) :: n

    if (n .eq. 0) then
      res = 1
    else
      res = n * SfManual_Factorial(n - 1)
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure recursive function SfManual_Binomial(n, k) result(res)

    integer(I32) :: res
    integer(I32), intent(in) :: n
    integer(I32), intent(in) :: k

    if (n .eq. k) then
      res = 1
    else if (k .eq. 0) then
      res = 1
    else if (n > k) then
      res = SfManual_Binomial(n - 1, k - 1) + SfManual_Binomial(n - 1, k)
    else if (n < k .and. n >= 0) then
      res = 0
    else
      error stop "error stop: Binomial undefined"
    end if

  end function

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  recursive function SfManual_Hermite(x, n) result(res)
    real(R64)                :: res
    real(R64), intent(in)    :: x
    integer(I32), intent(in) :: n

    if (n < 0) then
      res = 0.0_R64
    else if (n .eq. 0) then
      res = 1.0_R64
    else if (n .eq. 1) then
      res = 2.0_R64 * x
    else
      res = 2.0_R64 * x * SfManual_Hermite(x, n - 1) - &
            2.0_R64 * real(n - 1, R64) * SfManual_Hermite(x, n - 2)
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  recursive function SfManual_Laguerre(n, k, x) result(res)
    real(R64)                :: res
    integer(I32), intent(in) :: n
    integer(I32), intent(in) :: k
    real(R64), intent(in)    :: x

    if (n < 0) then
      res = 0.0_R64
    else if (n .eq. 0) then
      res = 1.0_R64
    else if (n .eq. 1) then
      res = 1.0_R64 + real(k, R64) - x
    else
      res = ((2.0_R64 * real(n - 1, R64) + real(k, R64) + 1.0_R64 - x) * SfManual_Laguerre(n - 1, k, x) - &
             (real(n - 1, R64) + real(k, R64)) * SfManual_Laguerre(n - 2, k, x)) / real(n, R64)
    end if

  end function

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  recursive function SfManual_Legendre(l, m, x) result(res)
    ! Computes the associated Legendre polynomial P_l^m(x)
    ! Uses standard recurrence relations.
    ! Note: m can be negative, P_l^{-m}(x) = (-1)^m * (l-m)! / (l+m)! * P_l^m(x)
    real(R64)                :: res
    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    real(R64), intent(in)    :: x
    real(R64)                :: somx2, fact
    integer(I32)             :: absM, i

    absM = abs(m)

    ! Error case
    if (l < 0 .or. absM > l) then
      res = 0.0_R64
      return
    end if

    ! Base case: P_m^m(x)
    if (l .eq. absM) then
      res = 1.0_R64
      if (absM > 0) then
        somx2 = sqrt(1.0_R64 - x**2)
        fact = 1.0_R64
        do i = 1, absM
          res = res * (-fact) * somx2
          fact = fact + 2.0_R64
        end do
      end if

      ! Base case: P_{m+1}^m(x)
    else if (l .eq. absM + 1) then
      res = x * (2.0_R64 * real(absM, R64) + 1.0_R64) * SfManual_Legendre(absM, absM, x)

      ! Recursive case using recurrence relation
    else
      res = ((2.0_R64 * real(l, R64) - 1.0_R64) * x * SfManual_Legendre(l - 1, absM, x) - &
             (real(l + absM, R64) - 1.0_R64) * SfManual_Legendre(l - 2, absM, x)) / &
            real(l - absM, R64)
    end if

    ! Handle negative m values
    if (m < 0) then
      fact = real(SfManual_Factorial(l - absM), R64) / real(SfManual_Factorial(l + absM), R64)
      if (mod(absM, 2) .eq. 1) then
        res = -res * fact
      else
        res = res * fact
      end if
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function SfManual_SpHarm(l, m, theta, phi) result(res)
    use M_Utils_Constants
    ! Computes the spherical harmonic Y_l^m(theta, phi)
    ! theta is the polar angle [0, pi], phi is the azimuthal angle [0, 2pi]
    complex(R64)             :: res
    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    real(R64), intent(in)    :: theta
    real(R64), intent(in)    :: phi
    real(R64)                :: factor, legendreVal, factRatio
    integer(I32)             :: absM

    absM = abs(m)

    factRatio = real(SfManual_Factorial(l - absM), R64) / real(SfManual_Factorial(l + absM), R64)

    factor = sqrt((2.0_R64 * real(l, R64) + 1.0_R64) / (4.0_R64 * PI) * factRatio)

    legendreVal = SfManual_Legendre(l, m, cos(theta))

    res = factor * legendreVal * exp(CMPLX(0.0_R64, real(m, R64) * phi, kind=R64))

  end function

end module
