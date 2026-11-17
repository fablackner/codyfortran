! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> SHTns (Spherical Harmonics Transform) thin wrappers and helpers.
!>
!> Provides convenience routines to configure SHTns contexts and perform
!> spectral transforms used by spherical-grid features.
module M_Utils_ShtnsLib

  use M_Utils_Types
  use, intrinsic :: iso_c_binding

  implicit none

  include 'shtns.f03'

  type :: T_ShtnsLib_Ctx
    type(C_PTR) :: handler = C_NULL_PTR
    type(shtns_info), pointer :: info => null()

    integer(I32) :: nlm = 0
    integer(I32) :: layout = 0
    integer(I32) :: norm = 0

    integer(I32) :: nlat = 0
    integer(I32) :: nphi = 0

    real(R64)    :: eps_polar = 0.0_R64
    logical      :: initializedQ = .false.

    integer(I32) :: lmax = 0
    integer(I32) :: mmax = 0
    integer(I32) :: mres = 0
  end type

  public :: T_ShtnsLib_Ctx

  public :: ShtnsLib_CreateCtx
  public :: ShtnsLib_DestroyCtx

  public :: ShtnsLib_GetGaussLegendreGrid
  public :: ShtnsLib_GetPhiGrid

  public :: ShtnsLib_DoSpHarmTransform
  public :: ShtnsLib_DoInverseSpHarmTransform

  public :: ShtnsLib_GetLmIndex
  public :: ShtnsLib_GetLM

  public :: ShtnsLib_GetNlm, ShtnsLib_GetLmax, ShtnsLib_GetMmax

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ShtnsLib_CreateCtx(ctx, lmax, mmax, mres, nlat, nphi, verboseFlag_)
    type(T_ShtnsLib_Ctx), intent(out) :: ctx
    integer(I32), intent(in) :: lmax, mmax, mres, nlat, nphi
    integer(I32), intent(in), optional :: verboseFlag_

    integer(I32) :: nthreads

    if (present(verboseFlag_)) then
      call shtns_verbose(verboseFlag_)
    else
      call shtns_verbose(0)
    end if

    nthreads = shtns_use_threads(0)

    ctx % lmax = lmax
    ctx % mmax = mmax
    ctx % mres = mres

    ctx % nlat = nlat
    ctx % nphi = nphi

    ctx % layout = SHT_PHI_CONTIGUOUS
    ctx % norm = SHT_ORTHONORMAL + SHT_REAL_NORM

    ctx % handler = shtns_create(lmax, mmax, mres, ctx % norm)
    if (.not. c_associated(ctx % handler)) error stop "ShtnsLib_CreateCtx: handler creation failed"

    call c_f_pointer(ctx % handler, ctx % info)

    ctx % nlm = ctx % info % nlm
    ctx % eps_polar = 1.0e-10_R64

    call shtns_set_grid(ctx % handler, SHT_GAUSS, ctx % eps_polar, nlat, nphi)

    ctx % nlat = ctx % info % nlat
    ctx % nphi = ctx % info % nphi

    if (present(verboseFlag_)) then
      if (verboseFlag_ > 0) call shtns_print_cfg(ctx % handler)
    end if

    ctx % initializedQ = .true.
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ShtnsLib_DestroyCtx(ctx)
    type(T_ShtnsLib_Ctx), intent(inout) :: ctx

    if (ctx % initializedQ .and. c_associated(ctx % handler)) then
      call shtns_destroy(ctx % handler)
    end if

    ctx % handler = c_null_ptr
    ctx % info => null()
    ctx % initializedQ = .false.
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ShtnsLib_GetGaussLegendreGrid(cos_theta, sin_theta, weights, ctx, theta_)
    real(R64), intent(out), allocatable :: cos_theta(:), sin_theta(:), weights(:)
    type(T_ShtnsLib_Ctx), intent(in) :: ctx
    real(R64), intent(out), allocatable, optional :: theta_(:)

    real(R64), pointer :: cos_ptr(:), sin_ptr(:)
    integer(I32) :: i

    if (.not. ctx % initializedQ) error stop "ShtnsLib_GetGaussLegendreGrid: ctx not initializedQ"

    allocate (cos_theta(ctx % nlat), sin_theta(ctx % nlat), weights(ctx % nlat))

    if (present(theta_)) allocate (theta_(ctx % nlat))

    call shtns_gauss_wts(ctx % handler, weights)

    call c_f_pointer(ctx % info % ct, cos_ptr, [ctx % nlat])
    call c_f_pointer(ctx % info % st, sin_ptr, [ctx % nlat])

    cos_theta = cos_ptr
    sin_theta = sin_ptr

    if (present(theta_)) then
      do i = 1, ctx % nlat
        theta_(i) = acos(cos_theta(i))
      end do
    end if
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ShtnsLib_GetPhiGrid(phi, ctx)
    use M_Utils_Constants

    real(R64), intent(out), allocatable :: phi(:)
    type(T_ShtnsLib_Ctx), intent(in) :: ctx

    integer(I32) :: i

    if (.not. ctx % initializedQ) error stop "ShtnsLib_GetPhiGrid: ctx not initializedQ"

    allocate (phi(ctx % nphi))

    do i = 1, ctx % nphi
      phi(i) = 2.0_R64 * PI * real(i - 1, R64) / real(ctx % nphi, R64)
    end do
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ShtnsLib_GetLmIndex(idx, l, m, ctx)
    integer(I32), intent(out) :: idx
    integer(I32), intent(in)  :: l, m
    type(T_ShtnsLib_Ctx), intent(in) :: ctx

    if (.not. ctx % initializedQ) error stop "ShtnsLib_GetLmIndex: ctx not initializedQ"

    idx = shtns_lmidx(ctx % handler, l, m)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ShtnsLib_GetLM(l, m, idx, ctx)
    integer(I32), intent(out) :: l, m
    integer(I32), intent(in)  :: idx
    type(T_ShtnsLib_Ctx), intent(in) :: ctx

    integer(I32) :: ll, mm

    if (.not. ctx % initializedQ) error stop "ShtnsLib_GetLM: ctx not initializedQ"

    do ll = 0, ctx % lmax
      do mm = -ll, ll, ctx % mres
        if (abs(mm) <= ctx % mmax) then
          if (shtns_lmidx(ctx % handler, ll, mm) .eq. idx) then
            l = ll
            m = mm
            return
          end if
        end if
      end do
    end do

    error stop "ShtnsLib_GetLM: invalid index"
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function ShtnsLib_GetNlm(ctx) result(res)
    integer(I32) :: res
    type(T_ShtnsLib_Ctx), intent(in) :: ctx

    if (.not. ctx % initializedQ) error stop "ShtnsLib_GetNlm: ctx not initializedQ"

    res = ctx % nlm
  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function ShtnsLib_GetLmax(ctx) result(res)
    integer(I32) :: res
    type(T_ShtnsLib_Ctx), intent(in) :: ctx

    if (.not. ctx % initializedQ) error stop "ShtnsLib_GetLmax: ctx not initializedQ"

    res = ctx % lmax
  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function ShtnsLib_GetMmax(ctx) result(res)
    type(T_ShtnsLib_Ctx), intent(in) :: ctx
    integer(I32) :: res

    if (.not. ctx % initializedQ) error stop "ShtnsLib_GetMmax: ctx not initializedQ"

    res = ctx % mmax
  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ShtnsLib_DoSpHarmTransform(sh_coeffs, spatial_field, ctx)
    complex(R64), intent(out) :: sh_coeffs(:)
    real(R64), intent(inout)  :: spatial_field(:, :)
    type(T_ShtnsLib_Ctx), intent(in) :: ctx

    if (.not. ctx % initializedQ) error stop "ShtnsLib_SphericalHarmTransform: ctx not initializedQ"

    if (size(spatial_field, 1) .ne. ctx % nphi .or. size(spatial_field, 2) .ne. ctx % nlat) &
      error stop "ShtnsLib_SphericalHarmTransform: spatial_field shape mismatch"

    if (size(sh_coeffs) .ne. ctx % nlm) error stop "ShtnsLib_SphericalHarmTransform: sh_coeffs size mismatch"

    call spat_to_SH(ctx % handler, spatial_field, sh_coeffs)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ShtnsLib_DoInverseSpHarmTransform(spatial_field, sh_coeffs, ctx)
    real(R64), intent(out)    :: spatial_field(:, :)
    complex(R64), intent(in)  :: sh_coeffs(:)
    type(T_ShtnsLib_Ctx), intent(in) :: ctx

    if (.not. ctx % initializedQ) error stop "ShtnsLib_InverseSphericalHarmTransform: ctx not initializedQ"

    if (size(spatial_field, 1) .ne. ctx % nphi .or. size(spatial_field, 2) .ne. ctx % nlat) &
      error stop "ShtnsLib_InverseSphericalHarmTransform: spatial_field shape mismatch"

    if (size(sh_coeffs) .ne. ctx % nlm) error stop "ShtnsLib_InverseSphericalHarmTransform: sh_coeffs size mismatch"

    call SH_to_spat(ctx % handler, sh_coeffs, spatial_field)
  end subroutine
end module
