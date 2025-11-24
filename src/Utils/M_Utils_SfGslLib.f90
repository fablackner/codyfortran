! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fortran interface to selected GSL special functions.
module M_Utils_SfGslLib
  use M_Utils_Types
  use, intrinsic :: iso_c_binding

  implicit none

  private :: GSL_SUCCESS, GSL_FAILURE
  integer(c_int), parameter :: GSL_SUCCESS = 0
  integer(c_int), parameter :: GSL_FAILURE = -1

  interface
    pure function GSL_SF_CHOOSE(n, m) bind(C, name="gsl_sf_choose") result(res)
      use, intrinsic :: iso_c_binding
      real(c_double) :: res
      integer(c_int), intent(in), value :: n
      integer(c_int), intent(in), value :: m
    end function

    pure function GSL_SF_FACT(n) bind(C, name="gsl_sf_fact") result(res)
      use, intrinsic :: iso_c_binding
      real(c_double) :: res
      integer(c_int), intent(in), value :: n
    end function

    pure function GSL_SF_HERMITE_FUNC(n, x) bind(C, name="gsl_sf_hermite_func") result(res)
      use, intrinsic :: iso_c_binding
      real(c_double) :: res
      integer(c_int), intent(in), value :: n
      real(c_double), intent(in), value :: x
    end function

    pure function GSL_SF_LAGUERRE_N(n, a, x) bind(C, name="gsl_sf_laguerre_n") result(res)
      use, intrinsic :: iso_c_binding
      real(c_double) :: res
      integer(c_int), intent(in), value :: n
      real(c_double), intent(in), value :: a
      real(c_double), intent(in), value :: x
    end function

    pure function GSL_SF_LEGENDRE_PLM(l, m, x) bind(C, name="gsl_sf_legendre_Plm") result(res)
      use, intrinsic :: iso_c_binding
      real(c_double) :: res
      integer(c_int), intent(in), value :: l
      integer(c_int), intent(in), value :: m
      real(c_double), intent(in), value :: x
    end function

    pure function GSL_SF_LEGENDRE_SPHPLM(l, m, x) bind(C, name="gsl_sf_legendre_sphPlm") result(res)
      use, intrinsic :: iso_c_binding
      real(c_double) :: res
      integer(c_int), intent(in), value :: l
      integer(c_int), intent(in), value :: m
      real(c_double), intent(in), value :: x
    end function

    pure function GSL_SF_COUPLING_3J(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc) &
      bind(C, name="gsl_sf_coupling_3j") result(res)
      use, intrinsic :: iso_c_binding
      real(c_double) :: res
      integer(c_int), intent(in), value :: two_ja, two_jb, two_jc
      integer(c_int), intent(in), value :: two_ma, two_mb, two_mc
    end function

  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function SfGslLib_Factorial(n) result(res)
    use, intrinsic :: iso_c_binding

    integer(I32) :: res
    integer(I32), intent(in) :: n

    integer(c_int) :: nGsl
    real(c_double) :: resGsl

    nGsl = n

    if (n >= 0) then
      resGsl = GSL_SF_FACT(nGsl)
      res = nint(resGsl, kind=I32)
    else
      error stop "error stop: Factorial undefined"
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function SfGslLib_Binomial(n, k) result(res)
    use, intrinsic :: iso_c_binding

    integer(I32) :: res
    integer(I32), intent(in) :: n
    integer(I32), intent(in) :: k

    integer(c_int) :: nGsl, kGsl
    real(c_double) :: resGsl

    nGsl = n
    kGsl = k

    if (n .eq. k) then
      res = 1
    else if (k .eq. 0) then
      res = 1
    else if (n > k) then
      resGsl = GSL_SF_CHOOSE(nGsl, kGsl)
      res = nint(resGsl, kind=I32)
    else if (n < k .and. n >= 0) then
      res = 0
    else
      error stop "error stop: Binomial undefined"
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function SfGslLib_Hermite(x, n) result(res)
    use M_Utils_Constants

    real(R64)                :: res
    real(R64), intent(in)    :: x
    integer(I32), intent(in) :: n

    integer(c_int) :: n_c
    real(c_double) :: x_c, res_c

    if (n < 0) then
      res = 0.0_R64
      return
    end if

    n_c = n
    x_c = x
    res_c = GSL_SF_HERMITE_FUNC(n_c, x_c)

    ! GSL's hermite_func returns a normalized version, multiply by normalization factor
    res = res_c * sqrt(2.0_R64**n * SfGslLib_Factorial(n) * sqrt(PI))

  end function

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function SfGslLib_Laguerre(n, k, x) result(res)
    real(R64)                :: res
    integer(I32), intent(in) :: n
    integer(I32), intent(in) :: k
    real(R64), intent(in)    :: x

    integer(c_int) :: n_c
    real(c_double) :: k_c, x_c, res_c

    if (n < 0) then
      res = 0.0_R64
      return
    end if

    n_c = n
    k_c = real(k, c_double)
    x_c = x

    res_c = GSL_SF_LAGUERRE_N(n_c, k_c, x_c)
    res = res_c

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function SfGslLib_Legendre(l, m, x) result(res)
    real(R64)                :: res
    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    real(R64), intent(in)    :: x

    integer(c_int) :: l_c, m_c
    real(c_double) :: x_c, res_c

    l_c = l
    m_c = m
    x_c = x

    res_c = GSL_SF_LEGENDRE_PLM(l_c, abs(m_c), x_c)

    ! Handle negative m values
    if (m < 0) then
      res = (-1)**abs(m) * (SfGslLib_Factorial(l - abs(m)) / SfGslLib_Factorial(l + abs(m))) * res_c
    else
      res = res_c
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function SfGslLib_SpHarm(l, m, theta, phi) result(res)
    complex(R64)             :: res
    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    real(R64), intent(in)    :: theta
    real(R64), intent(in)    :: phi

    integer(c_int) :: l_c, m_c
    real(c_double) :: cos_theta_c, res_c

    l_c = l
    m_c = m
    cos_theta_c = cos(theta)

    ! Get the normalized spherical harmonic (without the azimuthal part)
    res_c = GSL_SF_LEGENDRE_SPHPLM(l_c, abs(m_c), cos_theta_c)

    ! Handle negative m values and add azimuthal part
    if (m < 0) then
      res = (-1)**abs(m) * cmplx(res_c, 0.0_R64, kind=R64) * exp(cmplx(0.0_R64, real(m, R64) * phi, kind=R64))
    else
      res = cmplx(res_c, 0.0_R64, kind=R64) * exp(cmplx(0.0_R64, real(m, R64) * phi, kind=R64))
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function SfGslLib_Wigner3jGsl(j1, j2, j3, m1, m2, m3) result(res)
    ! Computes the Wigner 3j symbol using GSL:
    !
    ! / j1  j2  j3 \
    ! \ m1  m2  m3 /
    !
    ! GSL expects integer arguments that are twice the actual j and m values
    ! because angular momentum can be half-integer in quantum mechanics.
    ! For our purposes with spherical harmonics, j and m are integers.

    real(R64) :: res
    integer(I32), intent(in) :: j1, j2, j3, m1, m2, m3

    ! Call GSL's 3j symbol function
    res = GSL_SF_COUPLING_3J(2 * j1, 2 * j2, 2 * j3, 2 * m1, 2 * m2, 2 * m3)

  end function

end module
