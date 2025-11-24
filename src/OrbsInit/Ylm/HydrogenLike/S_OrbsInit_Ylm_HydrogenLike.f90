! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_OrbsInit_Ylm_HydrogenLike) S_OrbsInit_Ylm_HydrogenLike

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine OrbsInit_Ylm_HydrogenLike_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit_Ylm

    implicit none

    real(R64) :: charge

    call Say_Fabricate("orbsInit.ylm.hydrogenLike")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    charge = Json_Get("orbsInit.ylm.hydrogenLike.charge", 1.0_R64)
    OrbsInit_Ylm_HydrogenLike_charge = charge

    ! Read sequences for quantum numbers
    OrbsInit_Ylm_HydrogenLike_n = Json_Get("orbsInit.ylm.hydrogenLike.n", [1])
    OrbsInit_Ylm_HydrogenLike_l = Json_Get("orbsInit.ylm.hydrogenLike.l", [0])
    OrbsInit_Ylm_HydrogenLike_m = Json_Get("orbsInit.ylm.hydrogenLike.m", [0])

    OrbsInit_Ylm_InitFunction => InitFunction

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function InitFunction(r, l, m, index, bt_) result(res)
    use M_Utils_UnusedVariables
    use M_Utils_SfGslLib

    complex(R64)             :: res
    real(R64), intent(in)    :: r
    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    integer(I32), intent(in) :: index
    integer(I32), intent(in), optional :: bt_

    real(R64) :: charge, a0
    integer(I32) :: n
    real(R64)    :: rho

    if (.false.) call UnusedVariables_Mark(m, bt_)

    res = (0.0_R64, 0.0_R64)

    if (l .ne. OrbsInit_Ylm_HydrogenLike_l(index)) return
    if (m .ne. OrbsInit_Ylm_HydrogenLike_m(index)) return
    n = OrbsInit_Ylm_HydrogenLike_n(index)

    ! Get parameters
    charge = OrbsInit_Ylm_HydrogenLike_charge
    a0 = 1.0_R64    ! Bohr radius (in atomic units)

    ! Calculate rho parameter
    rho = 2.0_R64 * charge * r / (n * a0)

    ! Calculate radial part R_nl(r) for hydrogen-like atoms
    res = exp(-rho / 2.0_R64) * (rho**l) * &
          SfGslLib_Laguerre(n - l - 1, 2 * l + 1, rho)

  end function

end submodule
