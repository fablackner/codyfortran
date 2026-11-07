! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Exchange-correlation functionals of Libxc (C API, version 5 or later).
!>
!> A functional is an opaque Libxc handle (`xc_func_alloc`) inside
!> `T_XcLib_Functional`, so nothing here depends on the layout of Libxc's
!> structs. Evaluation is on point arrays in Libxc's layout:
!>
!> | nSpin | rho            | sigma                       | vrho      | vsigma                      |
!> |-------|----------------|-----------------------------|-----------|-----------------------------|
!> | 1     | rho(np)        | sigma(np) = \|∇ρ\|²          | vrho(np)  | vsigma(np)                  |
!> | 2     | rho(2, np): ↑ ↓ | sigma(3, np): ∇ρ↑·∇ρ↑, ∇ρ↑·∇ρ↓, ∇ρ↓·∇ρ↓ | vrho(2, np) | vsigma(3, np) |
!>
!> zk(np) is the energy per particle, E_xc = ∫ (ρ↑ + ρ↓) zk, and vrho, vsigma
!> are the partial derivatives of the energy density (ρ↑ + ρ↓) zk. The
!> evaluation routines only read the handle and are reentrant: they may be
!> called concurrently from OpenMP threads on disjoint chunks of points.
!>
!> The implementation submodule is chosen at build time: S_Utils_XcLib_Libxc
!> (CMake option CODY_USE_LIBXC=ON, the default; links libxc through
!> bind(C)) or S_Utils_XcLib_Stub (no Libxc: `XcLib_AvailableQ` is false and
!> every other call stops).
module M_Utils_XcLib
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! constants of xc.h
  !=============================================================================

  !> Spin channels of `XcLib_Init`
  integer(I32), parameter :: XCLIB_UNPOLARIZED = 1
  integer(I32), parameter :: XCLIB_POLARIZED = 2

  !> Kinds
  integer(I32), parameter :: XCLIB_EXCHANGE = 0
  integer(I32), parameter :: XCLIB_CORRELATION = 1
  integer(I32), parameter :: XCLIB_EXCHANGE_CORRELATION = 2
  integer(I32), parameter :: XCLIB_KINETIC = 3

  !> Families (the hybrid families carry a fraction of exact exchange)
  integer(I32), parameter :: XCLIB_FAMILY_LDA = 1
  integer(I32), parameter :: XCLIB_FAMILY_GGA = 2
  integer(I32), parameter :: XCLIB_FAMILY_MGGA = 4
  integer(I32), parameter :: XCLIB_FAMILY_LCA = 8
  integer(I32), parameter :: XCLIB_FAMILY_OEP = 16
  integer(I32), parameter :: XCLIB_FAMILY_HYB_GGA = 32
  integer(I32), parameter :: XCLIB_FAMILY_HYB_MGGA = 64
  integer(I32), parameter :: XCLIB_FAMILY_HYB_LDA = 128

  !> Flags (bits of `T_XcLib_Functional % flags`)
  integer(I32), parameter :: XCLIB_FLAGS_HAVE_EXC = 1
  integer(I32), parameter :: XCLIB_FLAGS_HAVE_VXC = 2
  integer(I32), parameter :: XCLIB_FLAGS_1D = 32
  integer(I32), parameter :: XCLIB_FLAGS_2D = 64
  integer(I32), parameter :: XCLIB_FLAGS_3D = 128
  integer(I32), parameter :: XCLIB_FLAGS_HYB_CAM = 256
  integer(I32), parameter :: XCLIB_FLAGS_HYB_CAMY = 512
  integer(I32), parameter :: XCLIB_FLAGS_VV10 = 1024
  integer(I32), parameter :: XCLIB_FLAGS_HYB_LC = 2048
  integer(I32), parameter :: XCLIB_FLAGS_HYB_LCY = 4096
  integer(I32), parameter :: XCLIB_FLAGS_NEEDS_LAPLACIAN = 32768
  integer(I32), parameter :: XCLIB_FLAGS_NEEDS_TAU = 65536

  !=============================================================================
  ! types
  !=============================================================================

  !> An initialized Libxc functional and its properties (`XcLib_Init`).
  type :: T_XcLib_Functional
    !> Libxc handle (xc_func_type *), null before `XcLib_Init` and after `XcLib_End`
    type(c_ptr) :: handle = c_null_ptr
    !> Libxc functional number, e.g. 101 for gga_x_pbe
    integer(I32) :: id = 0
    !> Spin channels: XCLIB_UNPOLARIZED or XCLIB_POLARIZED
    integer(I32) :: nSpin = 0
    !> Family (XCLIB_FAMILY_*)
    integer(I32) :: family = 0
    !> Kind (XCLIB_EXCHANGE, ...)
    integer(I32) :: kind = 0
    !> Flags (XCLIB_FLAGS_* bits)
    integer(I32) :: flags = 0
    !> Number of external parameters
    integer(I32) :: nExtParams = 0
    !> Descriptive name of Libxc, e.g. "Perdew, Burke & Ernzerhof"
    character(len=:), allocatable :: description
  end type

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> True if CodyFortran was built with Libxc (CODY_USE_LIBXC=ON).
    module function XcLib_AvailableQ() result(res)
      logical :: res
    end function

    !> Version string of the linked Libxc, e.g. "7.0.0".
    module function XcLib_Version() result(res)
      character(len=:), allocatable :: res
    end function

    !> Libxc number of a functional name such as "gga_x_pbe" or "XC_GGA_X_PBE"
    !> (case-insensitive), -1 if Libxc does not know the name.
    module function XcLib_FunctionalId(name) result(id)
      character(len=*), intent(in) :: name
      integer(I32) :: id
    end function

    !> Allocates and initializes functional `id` for nSpin channels and fills
    !> its properties; stops if Libxc cannot initialize it.
    module subroutine XcLib_Init(func, id, nSpin)
      type(T_XcLib_Functional), intent(out) :: func
      integer(I32), intent(in) :: id
      integer(I32), intent(in) :: nSpin
    end subroutine

    !> Releases the Libxc handle (no-op for an uninitialized functional).
    module subroutine XcLib_End(func)
      type(T_XcLib_Functional), intent(inout) :: func
    end subroutine

    !> Densities below threshold are treated as zero by the functional.
    module subroutine XcLib_SetDensityThreshold(func, threshold)
      type(T_XcLib_Functional), intent(inout) :: func
      real(R64), intent(in) :: threshold
    end subroutine

    !> Sets all external parameters (size(params) must be func % nExtParams).
    module subroutine XcLib_SetExtParams(func, params)
      type(T_XcLib_Functional), intent(inout) :: func
      real(R64), intent(in), contiguous :: params(:)
    end subroutine

    !> Default values of the external parameters, in Libxc's order.
    module subroutine XcLib_FillExtParamsDefault(params, func)
      real(R64), intent(out), allocatable :: params(:)
      type(T_XcLib_Functional), intent(in) :: func
    end subroutine

    !> Fraction of exact (Hartree-Fock) exchange of a global hybrid
    !> (xc_hyb_exx_coef); 0 for the semilocal families.
    module function XcLib_ExactExchangeFraction(func) result(res)
      type(T_XcLib_Functional), intent(in) :: func
      real(R64) :: res
    end function

    !> LDA energy per particle and potential of size(zk) points
    !> (xc_lda_exc_vxc; layout in the module header).
    module subroutine XcLib_EvaluateLda(zk, vrho, func, rho)
      real(R64), intent(out), contiguous :: zk(:)
      real(R64), intent(out), contiguous :: vrho(:)
      type(T_XcLib_Functional), intent(in) :: func
      real(R64), intent(in), contiguous :: rho(:)
    end subroutine

    !> GGA energy per particle and derivatives of size(zk) points
    !> (xc_gga_exc_vxc; layout in the module header).
    module subroutine XcLib_EvaluateGga(zk, vrho, vsigma, func, rho, sigma)
      real(R64), intent(out), contiguous :: zk(:)
      real(R64), intent(out), contiguous :: vrho(:)
      real(R64), intent(out), contiguous :: vsigma(:)
      type(T_XcLib_Functional), intent(in) :: func
      real(R64), intent(in), contiguous :: rho(:)
      real(R64), intent(in), contiguous :: sigma(:)
    end subroutine
  end interface

end module
