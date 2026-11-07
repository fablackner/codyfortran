! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Stub implementation of M_Utils_XcLib for builds without Libxc (CMake
!> option CODY_USE_LIBXC=OFF): no functional is available and every call
!> except `XcLib_AvailableQ`, `XcLib_Version` and `XcLib_End` stops.
submodule(M_Utils_XcLib) S_Utils_XcLib_Stub
  use M_Utils_UnusedVariables

  implicit none

  character(len=*), parameter :: NO_LIBXC = "XcLib: CodyFortran was built without Libxc (cmake -DCODY_USE_LIBXC=ON)"

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function XcLib_AvailableQ() result(res)
    logical :: res

    res = .false.

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function XcLib_Version() result(res)
    character(len=:), allocatable :: res

    res = "none"

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function XcLib_FunctionalId(name) result(id)
    character(len=*), intent(in) :: name
    integer(I32) :: id

    if (.false.) call UnusedVariables_Mark(name)
    id = -1
    error stop NO_LIBXC

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_Init(func, id, nSpin)
    type(T_XcLib_Functional), intent(out) :: func
    integer(I32), intent(in) :: id
    integer(I32), intent(in) :: nSpin

    if (.false.) call UnusedVariables_Mark(id, nSpin)
    func % handle = c_null_ptr
    error stop NO_LIBXC

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_End(func)
    type(T_XcLib_Functional), intent(inout) :: func

    func % handle = c_null_ptr

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_SetDensityThreshold(func, threshold)
    type(T_XcLib_Functional), intent(inout) :: func
    real(R64), intent(in) :: threshold

    if (.false.) call UnusedVariables_Mark(threshold)
    func % handle = c_null_ptr
    error stop NO_LIBXC

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_SetExtParams(func, params)
    type(T_XcLib_Functional), intent(inout) :: func
    real(R64), intent(in), contiguous :: params(:)

    if (.false.) call UnusedVariables_Mark(params)
    func % handle = c_null_ptr
    error stop NO_LIBXC

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_FillExtParamsDefault(params, func)
    real(R64), intent(out), allocatable :: params(:)
    type(T_XcLib_Functional), intent(in) :: func

    if (.false.) call UnusedVariables_Mark(func % id)
    allocate (params(0))
    error stop NO_LIBXC

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function XcLib_ExactExchangeFraction(func) result(res)
    type(T_XcLib_Functional), intent(in) :: func
    real(R64) :: res

    if (.false.) call UnusedVariables_Mark(func % id)
    res = 0.0_R64
    error stop NO_LIBXC

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_EvaluateLda(zk, vrho, func, rho)
    real(R64), intent(out), contiguous :: zk(:)
    real(R64), intent(out), contiguous :: vrho(:)
    type(T_XcLib_Functional), intent(in) :: func
    real(R64), intent(in), contiguous :: rho(:)

    if (.false.) call UnusedVariables_Mark(func % id, rho)
    zk(:) = 0.0_R64
    vrho(:) = 0.0_R64
    error stop NO_LIBXC

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_EvaluateGga(zk, vrho, vsigma, func, rho, sigma)
    real(R64), intent(out), contiguous :: zk(:)
    real(R64), intent(out), contiguous :: vrho(:)
    real(R64), intent(out), contiguous :: vsigma(:)
    type(T_XcLib_Functional), intent(in) :: func
    real(R64), intent(in), contiguous :: rho(:)
    real(R64), intent(in), contiguous :: sigma(:)

    if (.false.) call UnusedVariables_Mark(func % id, rho, sigma)
    zk(:) = 0.0_R64
    vrho(:) = 0.0_R64
    vsigma(:) = 0.0_R64
    error stop NO_LIBXC

  end subroutine

end submodule
