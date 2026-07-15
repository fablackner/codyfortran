! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Associated Legendre functions of the first and second kind.
!>
!> Provides P_l^m(x) on the cut (|x| <= 1, Ferrers functions with
!> Condon-Shortley phase) and off the cut (x > 1, Hobson-type with positive
!> (x^2-1)^(m/2) prefactor), as well as Q_l^m(x) for x > 1. These are the
!> building blocks of prolate-spheroidal (two-center) coordinate methods:
!> the eta expansion uses normalized on-cut P, while the xi Green's function
!> P_l^m(xi<) Q_l^m(xi>) uses the off-cut functions.
!>
!> Conventions
!> -----------
!> - On the cut:  P_l^m(x) = (-1)^m (1-x^2)^(m/2) d^m/dx^m P_l(x)
!> - Off the cut: P_l^m(x) = (x^2-1)^(m/2) d^m/dx^m P_l(x)
!>                Q_l^m(x) = (x^2-1)^(m/2) d^m/dx^m Q_l(x)
!> - Wronskian:   (x^2-1) [P_l^m(x) Q_l^m'(x) - P_l^m'(x) Q_l^m(x)]
!>                  = (-1)^(m+1) (l+m)!/(l-m)!   (x > 1)
!>
!> Q_l^m uses Miller's downward recurrence (Q is the minimal solution of the
!> degree recurrence for x > 1, so upward recurrence is unstable) normalized
!> to the closed-form Q_0^m(x).
module M_Utils_Legendre
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fill P(0:lmax) with the associated Legendre functions P_l^m(x) of fixed
!> order m for all degrees l = 0..lmax (entries with l < m are zero).
!> Works on the cut (|x| <= 1) and off the cut (x > 1); the degree recurrence
!> is identical, only the P_m^m seed differs between the two regions.
!> @param P     Output values P_l^m(x), indexed 0:lmax
!> @param lmax  Maximum degree
!> @param m     Order (m >= 0)
!> @param x     Argument
  pure subroutine Legendre_FillP(P, lmax, m, x)
    real(R64), intent(out) :: P(0:)
    integer(I32), intent(in) :: lmax
    integer(I32), intent(in) :: m
    real(R64), intent(in) :: x

    integer(I32) :: k, l
    real(R64) :: pmm

    P(0:lmax) = 0.0_R64
    if (m > lmax) return

    ! Seed P_m^m = (2m-1)!! * (1-x^2)^(m/2) with Condon-Shortley phase on the
    ! cut, respectively (2m-1)!! * (x^2-1)^(m/2) off the cut
    pmm = 1.0_R64
    do k = 1, m
      if (abs(x) <= 1.0_R64) then
        pmm = -pmm * (2 * k - 1) * sqrt(1.0_R64 - x * x)
      else
        pmm = pmm * (2 * k - 1) * sqrt(x * x - 1.0_R64)
      end if
    end do
    P(m) = pmm

    if (lmax .eq. m) return
    P(m + 1) = (2 * m + 1) * x * pmm

    ! Upward degree recurrence (stable for P in both regions)
    do l = m + 1, lmax - 1
      P(l + 1) = ((2 * l + 1) * x * P(l) - (l + m) * P(l - 1)) / (l - m + 1)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fill Q(0:lmax) with the associated Legendre functions of the second kind
!> Q_l^m(x) of fixed order m for all degrees l = 0..lmax, valid for x > 1.
!>
!> Uses Miller's algorithm: an unnormalized downward recurrence started well
!> above lmax (Q is the dominant solution in the downward direction), then
!> normalization to the closed-form Q_0^m(x). The start buffer is chosen from
!> the asymptotic decay rate exp(-l*acosh(x)) so that seed contamination is
!> below machine precision.
!> @param Q     Output values Q_l^m(x), indexed 0:lmax
!> @param lmax  Maximum degree
!> @param m     Order (m >= 0)
!> @param x     Argument (x > 1)
  pure subroutine Legendre_FillQ(Q, lmax, m, x)
    real(R64), intent(out) :: Q(0:)
    integer(I32), intent(in) :: lmax
    integer(I32), intent(in) :: m
    real(R64), intent(in) :: x

    integer(I32) :: l, lStart, buffer
    real(R64) :: alpha, scale
    real(R64) :: qNext, qCur, qPrev

    if (.not. (x > 1.0_R64 + 1.0e-12_R64)) then
      error stop "Legendre_FillQ: argument must satisfy x > 1 (Q diverges at x = 1)"
    end if

    ! Downward-recurrence contamination decays like exp(-2*alpha*buffer)
    alpha = log(x + sqrt(x * x - 1.0_R64))
    buffer = max(30, nint(45.0_R64 / alpha) + 1)
    lStart = lmax + buffer

    ! Buffer phase: pair descent with arbitrary seed, rescaled against overflow
    qNext = 0.0_R64
    qCur = 1.0e-30_R64
    do l = lStart, lmax + 1, -1
      qPrev = ((2 * l + 1) * x * qCur - (l - m + 1) * qNext) / (l + m)
      qNext = qCur
      qCur = qPrev
      if (abs(qCur) > 1.0e260_R64) then
        qCur = qCur * 1.0e-260_R64
        qNext = qNext * 1.0e-260_R64
      end if
    end do

    ! Stored phase: continue the same descent, keeping all values
    Q(lmax) = qCur
    if (lmax > 0) Q(lmax - 1) = ((2 * lmax + 1) * x * qCur - (lmax - m + 1) * qNext) / (lmax + m)
    do l = lmax - 1, 1, -1
      Q(l - 1) = ((2 * l + 1) * x * Q(l) - (l - m + 1) * Q(l + 1)) / (l + m)
    end do

    ! Normalize to the exact Q_0^m
    scale = Legendre_Q0(m, x) / Q(0)
    Q(0:lmax) = Q(0:lmax) * scale

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Closed-form Q_0^m(x) for x > 1:
!>   Q_0(x)   = (1/2) ln((x+1)/(x-1))
!>   Q_0^m(x) = (x^2-1)^(m/2) d^m/dx^m Q_0(x)
!>            = (x^2-1)^(m/2) * (-1)^(m-1) (m-1)!/2 * [1/(x+1)^m - 1/(x-1)^m]
  pure function Legendre_Q0(m, x) result(res)
    integer(I32), intent(in) :: m
    real(R64), intent(in) :: x
    real(R64) :: res

    integer(I32) :: k
    real(R64) :: factor

    if (m .eq. 0) then
      res = 0.5_R64 * log((x + 1.0_R64) / (x - 1.0_R64))
      return
    end if

    ! factor = (x^2-1)^(m/2) * (-1)^(m-1) (m-1)! / 2
    factor = 0.5_R64
    do k = 1, m - 1
      factor = -factor * k
    end do
    factor = factor * sqrt(x * x - 1.0_R64)**m

    res = factor * (1.0_R64 / (x + 1.0_R64)**m - 1.0_R64 / (x - 1.0_R64)**m)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Factorial ratio (l+m)!/(l-m)! computed as a product (no overflow for the
!> moderate l, m used in multipole expansions).
  pure function Legendre_FactorialRatio(l, m) result(res)
    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    real(R64) :: res

    integer(I32) :: k

    res = 1.0_R64
    do k = l - m + 1, l + m
      res = res * k
    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> L2 normalization factor on [-1, 1] for the on-cut P_l^m:
!>   integral of [Legendre_NormFactorP(l,m) * P_l^m(x)]^2 dx = 1
  pure function Legendre_NormFactorP(l, m) result(res)
    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    real(R64) :: res

    res = sqrt((2 * l + 1) / (2.0_R64 * Legendre_FactorialRatio(l, m)))

  end function

end module
