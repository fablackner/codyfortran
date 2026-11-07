! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Libxc implementation of M_Utils_XcLib (CMake option CODY_USE_LIBXC=ON):
!> bind(C) interfaces to the C API of libxc (no Fortran modules of Libxc, so
!> any Fortran compiler works with a C build of the library).
submodule(M_Utils_XcLib) S_Utils_XcLib_Libxc
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_size_t, c_char, c_null_char, c_associated, c_f_pointer

  implicit none

  interface
    function xc_func_alloc() bind(C, name="xc_func_alloc") result(p)
      import :: c_ptr
      type(c_ptr) :: p
    end function

    function xc_func_init(p, functional, nspin) bind(C, name="xc_func_init") result(ierr)
      import :: c_ptr, c_int
      type(c_ptr), value :: p
      integer(c_int), value :: functional
      integer(c_int), value :: nspin
      integer(c_int) :: ierr
    end function

    subroutine xc_func_end(p) bind(C, name="xc_func_end")
      import :: c_ptr
      type(c_ptr), value :: p
    end subroutine

    subroutine xc_func_free(p) bind(C, name="xc_func_free")
      import :: c_ptr
      type(c_ptr), value :: p
    end subroutine

    function xc_func_get_info(p) bind(C, name="xc_func_get_info") result(info)
      import :: c_ptr
      type(c_ptr), value :: p
      type(c_ptr) :: info
    end function

    function xc_func_info_get_kind(info) bind(C, name="xc_func_info_get_kind") result(res)
      import :: c_ptr, c_int
      type(c_ptr), value :: info
      integer(c_int) :: res
    end function

    function xc_func_info_get_family(info) bind(C, name="xc_func_info_get_family") result(res)
      import :: c_ptr, c_int
      type(c_ptr), value :: info
      integer(c_int) :: res
    end function

    function xc_func_info_get_flags(info) bind(C, name="xc_func_info_get_flags") result(res)
      import :: c_ptr, c_int
      type(c_ptr), value :: info
      integer(c_int) :: res
    end function

    function xc_func_info_get_name(info) bind(C, name="xc_func_info_get_name") result(res)
      import :: c_ptr
      type(c_ptr), value :: info
      type(c_ptr) :: res
    end function

    function xc_func_info_get_n_ext_params(info) bind(C, name="xc_func_info_get_n_ext_params") result(res)
      import :: c_ptr, c_int
      type(c_ptr), value :: info
      integer(c_int) :: res
    end function

    function xc_func_info_get_ext_params_default_value(info, number) &
      bind(C, name="xc_func_info_get_ext_params_default_value") result(res)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: info
      integer(c_int), value :: number
      real(c_double) :: res
    end function

    function xc_functional_get_number(name) bind(C, name="xc_functional_get_number") result(res)
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int) :: res
    end function

    function xc_version_string() bind(C, name="xc_version_string") result(res)
      import :: c_ptr
      type(c_ptr) :: res
    end function

    subroutine xc_func_set_dens_threshold(p, threshold) bind(C, name="xc_func_set_dens_threshold")
      import :: c_ptr, c_double
      type(c_ptr), value :: p
      real(c_double), value :: threshold
    end subroutine

    subroutine xc_func_set_ext_params(p, params) bind(C, name="xc_func_set_ext_params")
      import :: c_ptr, c_double
      type(c_ptr), value :: p
      real(c_double), intent(in) :: params(*)
    end subroutine

    function xc_hyb_exx_coef(p) bind(C, name="xc_hyb_exx_coef") result(res)
      import :: c_ptr, c_double
      type(c_ptr), value :: p
      real(c_double) :: res
    end function

    subroutine xc_lda_exc_vxc(p, np, rho, zk, vrho) bind(C, name="xc_lda_exc_vxc")
      import :: c_ptr, c_size_t, c_double
      type(c_ptr), value :: p
      integer(c_size_t), value :: np
      real(c_double), intent(in) :: rho(*)
      real(c_double), intent(out) :: zk(*)
      real(c_double), intent(out) :: vrho(*)
    end subroutine

    subroutine xc_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma) bind(C, name="xc_gga_exc_vxc")
      import :: c_ptr, c_size_t, c_double
      type(c_ptr), value :: p
      integer(c_size_t), value :: np
      real(c_double), intent(in) :: rho(*)
      real(c_double), intent(in) :: sigma(*)
      real(c_double), intent(out) :: zk(*)
      real(c_double), intent(out) :: vrho(*)
      real(c_double), intent(out) :: vsigma(*)
    end subroutine
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function XcLib_AvailableQ() result(res)
    logical :: res

    res = .true.

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function XcLib_Version() result(res)
    character(len=:), allocatable :: res

    res = CString(xc_version_string())

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function XcLib_FunctionalId(name) result(id)
    character(len=*), intent(in) :: name
    integer(I32) :: id

    id = int(xc_functional_get_number(trim(adjustl(name))//c_null_char), kind=I32)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_Init(func, id, nSpin)
    type(T_XcLib_Functional), intent(out) :: func
    integer(I32), intent(in) :: id
    integer(I32), intent(in) :: nSpin

    type(c_ptr) :: info
    character(len=16) :: idString

    if (nSpin .ne. XCLIB_UNPOLARIZED .and. nSpin .ne. XCLIB_POLARIZED) error stop "XcLib_Init: nSpin must be 1 or 2"

    write (idString, '(I0)') id
    func % handle = xc_func_alloc()
    if (.not. c_associated(func % handle)) error stop "XcLib_Init: xc_func_alloc failed"
    if (xc_func_init(func % handle, int(id, kind=c_int), int(nSpin, kind=c_int)) .ne. 0) then
      call xc_func_free(func % handle)
      func % handle = c_null_ptr
      error stop "XcLib_Init: Libxc cannot initialize functional "//trim(idString)
    end if

    info = xc_func_get_info(func % handle)
    func % id = id
    func % nSpin = nSpin
    func % family = int(xc_func_info_get_family(info), kind=I32)
    func % kind = int(xc_func_info_get_kind(info), kind=I32)
    func % flags = int(xc_func_info_get_flags(info), kind=I32)
    func % nExtParams = int(xc_func_info_get_n_ext_params(info), kind=I32)
    func % description = CString(xc_func_info_get_name(info))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_End(func)
    type(T_XcLib_Functional), intent(inout) :: func

    if (.not. c_associated(func % handle)) return
    call xc_func_end(func % handle)
    call xc_func_free(func % handle)
    func % handle = c_null_ptr

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_SetDensityThreshold(func, threshold)
    type(T_XcLib_Functional), intent(inout) :: func
    real(R64), intent(in) :: threshold

    call CheckInitialized(func, "XcLib_SetDensityThreshold")
    call xc_func_set_dens_threshold(func % handle, real(threshold, kind=c_double))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_SetExtParams(func, params)
    type(T_XcLib_Functional), intent(inout) :: func
    real(R64), intent(in), contiguous :: params(:)

    call CheckInitialized(func, "XcLib_SetExtParams")
    if (size(params) .ne. func % nExtParams) error stop "XcLib_SetExtParams: wrong number of external parameters"
    if (func % nExtParams .eq. 0) return
    call xc_func_set_ext_params(func % handle, real(params, kind=c_double))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_FillExtParamsDefault(params, func)
    real(R64), intent(out), allocatable :: params(:)
    type(T_XcLib_Functional), intent(in) :: func

    type(c_ptr) :: info
    integer(I32) :: i

    call CheckInitialized(func, "XcLib_FillExtParamsDefault")
    allocate (params(func % nExtParams))
    info = xc_func_get_info(func % handle)
    do i = 1, func % nExtParams
      params(i) = real(xc_func_info_get_ext_params_default_value(info, int(i - 1, kind=c_int)), kind=R64)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function XcLib_ExactExchangeFraction(func) result(res)
    type(T_XcLib_Functional), intent(in) :: func
    real(R64) :: res

    call CheckInitialized(func, "XcLib_ExactExchangeFraction")
    res = 0.0_R64
    if (func % family .eq. XCLIB_FAMILY_HYB_LDA .or. func % family .eq. XCLIB_FAMILY_HYB_GGA .or. &
        func % family .eq. XCLIB_FAMILY_HYB_MGGA) res = real(xc_hyb_exx_coef(func % handle), kind=R64)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_EvaluateLda(zk, vrho, func, rho)
    real(R64), intent(out), contiguous :: zk(:)
    real(R64), intent(out), contiguous :: vrho(:)
    type(T_XcLib_Functional), intent(in) :: func
    real(R64), intent(in), contiguous :: rho(:)

    integer(I64) :: np

    np = size(zk, kind=I64)
    if (size(rho, kind=I64) .ne. func % nSpin * np .or. size(vrho, kind=I64) .ne. func % nSpin * np) &
      error stop "XcLib_EvaluateLda: array sizes do not match the number of points and spin channels"
    if (np .eq. 0) return
    call xc_lda_exc_vxc(func % handle, int(np, kind=c_size_t), rho, zk, vrho)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine XcLib_EvaluateGga(zk, vrho, vsigma, func, rho, sigma)
    real(R64), intent(out), contiguous :: zk(:)
    real(R64), intent(out), contiguous :: vrho(:)
    real(R64), intent(out), contiguous :: vsigma(:)
    type(T_XcLib_Functional), intent(in) :: func
    real(R64), intent(in), contiguous :: rho(:)
    real(R64), intent(in), contiguous :: sigma(:)

    integer(I64) :: np, nSigma

    np = size(zk, kind=I64)
    nSigma = 1
    if (func % nSpin .eq. XCLIB_POLARIZED) nSigma = 3
    if (size(rho, kind=I64) .ne. func % nSpin * np .or. size(vrho, kind=I64) .ne. func % nSpin * np .or. &
        size(sigma, kind=I64) .ne. nSigma * np .or. size(vsigma, kind=I64) .ne. nSigma * np) &
      error stop "XcLib_EvaluateGga: array sizes do not match the number of points and spin channels"
    if (np .eq. 0) return
    call xc_gga_exc_vxc(func % handle, int(np, kind=c_size_t), rho, sigma, zk, vrho, vsigma)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Fortran copy of a NUL-terminated C string (empty for a null pointer).
  function CString(cPtr) result(res)
    type(c_ptr), intent(in) :: cPtr
    character(len=:), allocatable :: res

    character(kind=c_char), pointer :: chars(:)
    integer(I32) :: n, i

    if (.not. c_associated(cPtr)) then
      res = ""
      return
    end if

    ! Libxc's names are short; 4096 bounds the scan of a corrupt pointer
    call c_f_pointer(cPtr, chars, [4096])
    n = 0
    do while (n < 4096)
      if (chars(n + 1) .eq. c_null_char) exit
      n = n + 1
    end do

    allocate (character(len=n) :: res)
    do i = 1, n
      res(i:i) = chars(i)
    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine CheckInitialized(func, caller)
    type(T_XcLib_Functional), intent(in) :: func
    character(len=*), intent(in) :: caller

    if (.not. c_associated(func % handle)) error stop caller//": functional not initialized (XcLib_Init)"

  end subroutine

end submodule
