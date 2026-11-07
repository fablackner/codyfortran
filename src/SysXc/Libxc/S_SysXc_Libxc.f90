! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Libxc functional of SysXc: fabrication, checks, chunked evaluation.
submodule(M_SysXc_Libxc) S_SysXc_Libxc
  use M_Utils_XcLib

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Functionals, (nFunctionals, nSpin): column 1 unpolarized, column 2 polarized
  type(T_XcLib_Functional), allocatable :: funcs(:, :)

  !> Weights of the functionals
  real(R64), allocatable :: coefficients(:)

  !> Whether functional i is evaluated with the GGA routine
  logical, allocatable :: ggaQ(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysXc_Libxc_Fabricate(path)
    use M_Utils_Json
    use M_Utils_Say
    use M_SysXc

    character(len=*), intent(in) :: path

    character(len=:), allocatable :: names(:), name
    character(len=0) :: noNames(0)
    real(R64), allocatable :: params(:)
    real(R64) :: densityThreshold, exxSum
    integer(I32) :: nF, i, s, id
    logical :: notGivenQ

    call Say_Fabricate(path)

    if (.not. XcLib_AvailableQ()) error stop "sysXc.*.libxc: CodyFortran was built without Libxc (cmake -DCODY_USE_LIBXC=ON)"
    if (.not. Json_GetExistence(path//".functionals")) &
      error stop "sysXc.*.libxc.functionals is missing (a list of Libxc names, e.g. [""lda_x"", ""lda_c_pw""])"

    !------------------------------------
    ! read configuration
    !------------------------------------

    names = Json_Get("functionals", noNames, path_=path)
    nF = size(names)
    coefficients = Json_Get("coefficients", [(1.0_R64, i=1, nF)], path_=path)
    if (size(coefficients) .ne. nF) error stop "sysXc.*.libxc.coefficients needs one weight per functional"
    densityThreshold = Json_Get("densityThreshold", 0.0_R64, path_=path)

    !------------------------------------
    ! initialize and check every functional, unpolarized and polarized
    !------------------------------------

    if (allocated(funcs)) then
      do s = 1, size(funcs, 2)
        do i = 1, size(funcs, 1)
          call XcLib_End(funcs(i, s))
        end do
      end do
      deallocate (funcs)
    end if
    if (allocated(ggaQ)) deallocate (ggaQ)
    allocate (funcs(nF, 2), ggaQ(nF))

    write (*, '(1X, 2A)') "Libxc version ", XcLib_Version()
    exxSum = 0.0_R64
    do i = 1, nF
      name = trim(adjustl(names(i)))
      id = XcLib_FunctionalId(name)
      if (id < 0) error stop "sysXc.*.libxc.functionals: unknown Libxc functional "//name

      do s = 1, 2
        call XcLib_Init(funcs(i, s), id, s)
        if (densityThreshold > 0.0_R64) call XcLib_SetDensityThreshold(funcs(i, s), densityThreshold)
        if (Json_GetExistence(path//".extParams."//name)) then
          call XcLib_FillExtParamsDefault(params, funcs(i, s))
          params = Json_Get(name, params, path_=path//".extParams")
          if (size(params) .ne. funcs(i, s) % nExtParams) &
            error stop "sysXc.*.libxc.extParams."//name//": wrong number of external parameters"
          call XcLib_SetExtParams(funcs(i, s), params)
        end if
      end do

      call CheckFunctional(funcs(i, 1), name, SysXc_dimension)
      ggaQ(i) = funcs(i, 1) % family .eq. XCLIB_FAMILY_GGA .or. funcs(i, 1) % family .eq. XCLIB_FAMILY_HYB_GGA
      exxSum = exxSum + coefficients(i) * XcLib_ExactExchangeFraction(funcs(i, 1))

      write (*, '(1X, A, I0, 5A, F8.4, A, F8.4, A, L1)') "  functional ", id, ": ", name, " (", &
        funcs(i, 1) % description, "), coefficient ", coefficients(i), ", exact exchange ", &
        XcLib_ExactExchangeFraction(funcs(i, 1)), ", gradients ", ggaQ(i)
    end do

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysXc_gradientQ = any(ggaQ)
    ! Json_Get with alternativeQ_ leaves the result undefined if the key is absent
    SysXc_exactExchangeFraction = Json_Get("exactExchangeFraction", exxSum, alternativeQ_=notGivenQ, path_=path)
    if (notGivenQ) SysXc_exactExchangeFraction = exxSum
    write (*, '(1X, A, F8.4)') "  exact exchange fraction ", SysXc_exactExchangeFraction

    SysXc_EvaluateFunctional => EvaluateFunctional

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Stops for functionals that SysXc cannot evaluate.
  subroutine CheckFunctional(func, name, dimension)
    type(T_XcLib_Functional), intent(in) :: func
    character(len=*), intent(in) :: name
    integer(I32), intent(in) :: dimension

    integer(I32) :: dimensionFlag

    select case (func % family)
    case (XCLIB_FAMILY_LDA, XCLIB_FAMILY_GGA, XCLIB_FAMILY_HYB_LDA, XCLIB_FAMILY_HYB_GGA)
    case (XCLIB_FAMILY_MGGA, XCLIB_FAMILY_HYB_MGGA)
      error stop "sysXc.*.libxc: meta-GGA functionals are not supported ("//name//")"
    case default
      error stop "sysXc.*.libxc: unsupported functional family ("//name//")"
    end select

    if (func % kind .eq. XCLIB_KINETIC) error stop "sysXc.*.libxc: kinetic-energy functionals are not supported ("//name//")"
    if (iand(func % flags, XCLIB_FLAGS_HAVE_EXC) .eq. 0 .or. iand(func % flags, XCLIB_FLAGS_HAVE_VXC) .eq. 0) &
      error stop "sysXc.*.libxc: the functional has no energy or potential ("//name//")"
    if (iand(func % flags, ior(ior(XCLIB_FLAGS_HYB_CAM, XCLIB_FLAGS_HYB_CAMY), &
                               ior(XCLIB_FLAGS_HYB_LC, XCLIB_FLAGS_HYB_LCY))) .ne. 0) &
      error stop "sysXc.*.libxc: range-separated hybrids are not supported (only global hybrids) ("//name//")"
    if (iand(func % flags, XCLIB_FLAGS_VV10) .ne. 0) &
      error stop "sysXc.*.libxc: VV10 nonlocal correlation is not supported ("//name//")"
    if (iand(func % flags, ior(XCLIB_FLAGS_NEEDS_TAU, XCLIB_FLAGS_NEEDS_LAPLACIAN)) .ne. 0) &
      error stop "sysXc.*.libxc: functionals of tau or the Laplacian are not supported ("//name//")"

    select case (dimension)
    case (1)
      dimensionFlag = XCLIB_FLAGS_1D
    case (2)
      dimensionFlag = XCLIB_FLAGS_2D
    case default
      dimensionFlag = XCLIB_FLAGS_3D
    end select
    if (iand(func % flags, dimensionFlag) .eq. 0) &
      error stop "sysXc.*.libxc: the functional is for another spatial dimension than the grid ("//name//")"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Σ_i c_i of the energy densities (ρ zk) and derivatives of the
  !> functionals, chunk by chunk in Libxc's interleaved layout.
  subroutine EvaluateFunctional(energyDensity, vRho, vSigma, rho, sigma)
    !$ use omp_lib, only: omp_in_parallel

    real(R64), intent(out), contiguous :: energyDensity(:)
    real(R64), intent(out), contiguous :: vRho(:, :)
    real(R64), intent(out), contiguous :: vSigma(:, :)
    real(R64), intent(in), contiguous :: rho(:, :)
    real(R64), intent(in), contiguous :: sigma(:, :)

    real(R64), allocatable :: rhoC(:), sigmaC(:), zk(:), vRhoC(:), vSigmaC(:)
    integer(I32) :: nP, nSpin, nSigma, nChunks, iChunk, p0, nInChunk, p, j, s, iF
    logical :: anyGgaQ

    nP = size(rho, 1)
    nSpin = size(rho, 2)
    nSigma = size(sigma, 2)
    anyGgaQ = any(ggaQ)

    if (anyGgaQ .and. nSigma .ne. 2 * nSpin - 1) error stop "SysXc_EvaluateFunctional (libxc): sigma has the wrong shape"
    if (size(vSigma, 2) .ne. nSigma) error stop "SysXc_EvaluateFunctional (libxc): vSigma has the wrong shape"

    energyDensity(:) = 0.0_R64
    vRho(:, :) = 0.0_R64
    vSigma(:, :) = 0.0_R64
    if (size(funcs, 1) .eq. 0) return

    nChunks = (nP + SysXc_Libxc_CHUNK_SIZE - 1) / SysXc_Libxc_CHUNK_SIZE

    !$omp parallel default(shared) private(rhoC, sigmaC, zk, vRhoC, vSigmaC, iChunk, p0, nInChunk, p, j, s, iF) &
    !$omp   if (.not. omp_in_parallel() .and. nChunks > 1)
    allocate (rhoC(nSpin * SysXc_Libxc_CHUNK_SIZE), zk(SysXc_Libxc_CHUNK_SIZE), vRhoC(nSpin * SysXc_Libxc_CHUNK_SIZE))
    allocate (sigmaC(max(nSigma, 1) * SysXc_Libxc_CHUNK_SIZE), vSigmaC(max(nSigma, 1) * SysXc_Libxc_CHUNK_SIZE))
    !$omp do schedule(dynamic)
    do iChunk = 1, nChunks
      p0 = (iChunk - 1) * SysXc_Libxc_CHUNK_SIZE
      nInChunk = min(SysXc_Libxc_CHUNK_SIZE, nP - p0)

      ! Interleave the spin channels (point fastest in our layout, spin in Libxc's)
      do j = 1, nInChunk
        p = p0 + j
        do s = 1, nSpin
          rhoC(nSpin * (j - 1) + s) = rho(p, s)
        end do
        do s = 1, nSigma
          sigmaC(nSigma * (j - 1) + s) = sigma(p, s)
        end do
      end do

      do iF = 1, size(funcs, 1)
        if (ggaQ(iF)) then
          call XcLib_EvaluateGga(zk(1:nInChunk), vRhoC(1:nSpin * nInChunk), vSigmaC(1:nSigma * nInChunk), funcs(iF, nSpin), &
                                 rhoC(1:nSpin * nInChunk), sigmaC(1:nSigma * nInChunk))
        else
          call XcLib_EvaluateLda(zk(1:nInChunk), vRhoC(1:nSpin * nInChunk), funcs(iF, nSpin), rhoC(1:nSpin * nInChunk))
        end if

        do j = 1, nInChunk
          p = p0 + j
          energyDensity(p) = energyDensity(p) + coefficients(iF) * sum(rho(p, :)) * zk(j)
          do s = 1, nSpin
            vRho(p, s) = vRho(p, s) + coefficients(iF) * vRhoC(nSpin * (j - 1) + s)
          end do
          if (ggaQ(iF)) then
            do s = 1, nSigma
              vSigma(p, s) = vSigma(p, s) + coefficients(iF) * vSigmaC(nSigma * (j - 1) + s)
            end do
          end if
        end do
      end do
    end do
    !$omp end do
    deallocate (rhoC, zk, vRhoC, sigmaC, vSigmaC)
    !$omp end parallel

  end subroutine

end submodule
