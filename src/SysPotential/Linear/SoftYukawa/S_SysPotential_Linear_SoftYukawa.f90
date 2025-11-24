! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Linear_SoftYukawa) S_SysPotential_Linear_SoftYukawa

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Linear_SoftYukawa_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential_Linear_SoftYukawa_StdImpl

    implicit none

    real(R64), allocatable :: position(:)
    real(R64), allocatable :: charge(:)
    real(R64), allocatable :: softening1(:)
    real(R64), allocatable :: softening2(:)
    real(R64), allocatable :: dampening(:)

    call Say_Fabricate("sysPotential.linear.softYukawa")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    position = Json_Get("sysPotential.linear.softYukawa.position", [0.0_R64])
    charge = Json_Get("sysPotential.linear.softYukawa.charge", [1.0_R64])
    softening1 = Json_Get("sysPotential.linear.softYukawa.softening1", [1.0_R64])
    softening2 = Json_Get("sysPotential.linear.softYukawa.softening2", [0.0_R64])
    dampening = Json_Get("sysPotential.linear.softYukawa.dampening", [0.0_R64])

    SysPotential_Linear_SoftYukawa_position = position
    SysPotential_Linear_SoftYukawa_charge = charge
    SysPotential_Linear_SoftYukawa_softening1 = softening1
    SysPotential_Linear_SoftYukawa_softening2 = softening2
    SysPotential_Linear_SoftYukawa_dampening = dampening

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.linear.softYukawa.stdImpl")) then
      call SysPotential_Linear_SoftYukawa_StdImpl_Fabricate

    else
      error stop "sysPotential.linear.softYukawa is missing one of: stdImpl"
    end if

  end subroutine

end submodule
