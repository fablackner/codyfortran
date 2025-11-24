! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Linear_Harmonic) S_SysPotential_Linear_Harmonic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Linear_Harmonic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential_Linear_Harmonic_StdImpl

    implicit none

    real(R64), allocatable :: position(:)
    real(R64), allocatable :: omega(:)

    call Say_Fabricate("sysPotential.linear.harmonic")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    position = Json_Get("sysPotential.linear.harmonic.position", [0.0_R64])
    omega = Json_Get("sysPotential.linear.harmonic.omega", [1.0_R64])
    SysPotential_Linear_Harmonic_position = position
    SysPotential_Linear_Harmonic_omega = omega

    ! Make sure all arrays have the same size
    if (size(position) .ne. size(omega)) then
      error stop "sysPotential.linear.harmonic: position and omega arrays must have the same length"
    end if

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.linear.harmonic.stdImpl")) then
      call SysPotential_Linear_Harmonic_StdImpl_Fabricate

    else
      error stop "sysPotential.linear.harmonic is missing one of: stdImpl"
    end if

  end subroutine

end submodule
