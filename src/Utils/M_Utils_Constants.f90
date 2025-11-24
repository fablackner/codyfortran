! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Centralized mathematical and physical constants with high precision.
module M_Utils_Constants
  use M_Utils_Types

  implicit none

  ! real(R64), parameter    :: SQRT2 = 1.41421356237309504880168872420969807856967_R64
  ! real(R64), parameter    :: EULER = 0.5772156649015328606065120900824024310422_R64
  ! real(R64), parameter    :: PI = 3.141592653589793238462643383279502884197_R64
  ! real(R64), parameter    :: PIO2 = 1.57079632679489661923132169163975144209858_R64
  ! real(R64), parameter    :: TWOPI = 6.283185307179586476925286766559005768394_R64
  real(R64), parameter :: SQRT2 = 1.41421356237309505_R64
  real(R64), parameter :: EULER = 0.57721566490153286_R64
  real(R64), parameter :: PI = 3.14159265358979324_R64
  real(R64), parameter :: PIO2 = 1.57079632679489662_R64
  real(R64), parameter :: TWOPI = 6.28318530717958648_R64

  complex(R64), parameter :: IU = (0.0_R64, 1.0_R64)
  complex(R64), parameter :: CONE = (1.0_R64, 0.0_R64)
  complex(R64), parameter :: CZERO = (0.0_R64, 0.0_R64)

  real(R64), parameter :: AUas = 1.0_R64 / 24.188843265_R64     ! attosecond in a.u.
  real(R64), parameter :: AUwcm2toel2 = 1.0_R64 / 3.50944481e16_R64 ! W/cm^2 in a.u. for electric field squared
  real(R64), parameter :: AUwcm2 = 1.55366146e-16_R64          ! W/cm^2 in a.u
  real(R64), parameter :: AUm = 1.0_R64 / 5.291772098e-11_R64    ! m in a.u.
  real(R64), parameter :: AUcm = 1.0_R64 / 5.291772098e-9_R64     ! cm in a.u.
  real(R64), parameter :: AUnm = 1.0_R64 / 5.291772098e-2_R64     ! nm in a.u.
  real(R64), parameter :: AUc = 137.03599911_R64          ! c (speed of light) in a.u. .eq. 1/alpha
  real(R64), parameter :: AUeV = 1.0_R64 / 27.2113845_R64       ! eV in a.u.

  !Frequently used physical constants
  real(R64), parameter :: HBAR = 1.054571628e-34_R64 ! in SI
  real(R64), parameter :: AMU = 1.660538782e-27_R64 ! in kg
  real(R64), parameter :: AEU = 4.35974381e-18_R64 ! in J

end module
