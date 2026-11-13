!> ctest-requires: libxc
!> Libxc binding (M_Utils_XcLib): names and numbers, Slater exchange against
!> its closed form (unpolarized, polarized, fully polarized), the PBE
!> exchange at zero gradient (= Slater), v_sigma of PBE against a finite
!> difference of the energy density, hybrid fractions (PBE0, B3LYP) and the
!> external parameters of the 1D soft-Coulomb exchange.
program T_Utils_XcLib
  use M_Utils_Types
  use M_Utils_Constants
  use M_Utils_XcLib
  use testdrive, only: check, error_type
  implicit none

  type(error_type), allocatable :: error
  type(T_XcLib_Functional) :: funcU, funcP
  real(R64), allocatable :: params(:)
  real(R64) :: rho(3), rho2(6), sigma(3), sigma3(9), zk(3), vrho(3), vrho2(6), vsigma(3), vsigma3(9)
  real(R64) :: zkP(3), vrhoP(3), vsigmaP(3), zkM(3), vrhoM(3), vsigmaM(3)
  real(R64) :: cx, maxErr, h
  integer(I32) :: i

  write (*, *) "Libxc version ", XcLib_Version()
  call check(error, XcLib_AvailableQ())
  if (allocated(error)) error stop "T_Utils_XcLib: Libxc not available"

  !------------------------------------
  ! names
  !------------------------------------

  call check(error, XcLib_FunctionalId("lda_x"), 1)
  if (allocated(error)) error stop "T_Utils_XcLib: lda_x id"
  call check(error, XcLib_FunctionalId("XC_GGA_X_PBE"), 101)
  if (allocated(error)) error stop "T_Utils_XcLib: XC_GGA_X_PBE id"
  call check(error, XcLib_FunctionalId("no_such_functional"), -1)
  if (allocated(error)) error stop "T_Utils_XcLib: unknown name"

  !------------------------------------
  ! Slater exchange: zk = -(3/4)(3/π)^(1/3) ρ^(1/3), vrho = -(3/π)^(1/3) ρ^(1/3)
  !------------------------------------

  rho = [1.0e-3_R64, 0.3_R64, 7.0_R64]
  cx = (3.0_R64 / PI)**(1.0_R64 / 3.0_R64)

  call XcLib_Init(funcU, XcLib_FunctionalId("lda_x"), XCLIB_UNPOLARIZED)
  call check(error, funcU % family, XCLIB_FAMILY_LDA)
  if (allocated(error)) error stop "T_Utils_XcLib: lda_x family"
  call check(error, iand(funcU % flags, XCLIB_FLAGS_3D) .ne. 0)
  if (allocated(error)) error stop "T_Utils_XcLib: lda_x is 3D"
  call XcLib_EvaluateLda(zk, vrho, funcU, rho)
  maxErr = max(maxval(abs(zk + 0.75_R64 * cx * rho**(1.0_R64 / 3.0_R64)) / abs(zk)), &
               maxval(abs(vrho + cx * rho**(1.0_R64 / 3.0_R64)) / abs(vrho)))
  write (*, *) "Slater exchange, unpolarized, relative error ", maxErr
  call check(error, maxErr < 1.0e-13_R64)
  if (allocated(error)) error stop "T_Utils_XcLib: Slater exchange (unpolarized)"

  ! Polarized with equal spins: the same energy per particle, each spin the potential
  call XcLib_Init(funcP, XcLib_FunctionalId("lda_x"), XCLIB_POLARIZED)
  do i = 1, 3
    rho2(2 * i - 1:2 * i) = 0.5_R64 * rho(i)
  end do
  call XcLib_EvaluateLda(zkP, vrho2, funcP, rho2)
  maxErr = max(maxval(abs(zkP - zk)/abs(zk)), maxval(abs(vrho2(1::2) - vrho)/abs(vrho)), &
               maxval(abs(vrho2(2::2) - vrho) / abs(vrho)))
  write (*, *) "Slater exchange, polarized (equal spins) vs unpolarized, relative difference ", maxErr
  call check(error, maxErr < 1.0e-13_R64)
  if (allocated(error)) error stop "T_Utils_XcLib: Slater exchange (polarized, equal spins)"

  ! Fully polarized: zk(ρ, ζ = 1) = 2^(1/3) zk(ρ, ζ = 0)
  do i = 1, 3
    rho2(2 * i - 1) = rho(i)
    rho2(2 * i) = 0.0_R64
  end do
  call XcLib_EvaluateLda(zkP, vrho2, funcP, rho2)
  maxErr = maxval(abs(zkP - 2.0_R64**(1.0_R64 / 3.0_R64) * zk) / abs(zk))
  write (*, *) "Slater exchange, fully polarized, relative error ", maxErr
  ! Libxc clips zeta at 1 - zeta_threshold (relative effect ~1e-12)
  call check(error, maxErr < 1.0e-10_R64)
  if (allocated(error)) error stop "T_Utils_XcLib: Slater exchange (fully polarized)"
  call XcLib_End(funcP)

  !------------------------------------
  ! PBE exchange: F_x(s = 0) = 1 (Slater), and v_sigma = d(ρ zk)/dσ
  !------------------------------------

  call XcLib_End(funcU)
  call XcLib_Init(funcU, XcLib_FunctionalId("gga_x_pbe"), XCLIB_UNPOLARIZED)
  call check(error, funcU % family, XCLIB_FAMILY_GGA)
  if (allocated(error)) error stop "T_Utils_XcLib: gga_x_pbe family"
  sigma = 0.0_R64
  call XcLib_EvaluateGga(zkP, vrhoP, vsigmaP, funcU, rho, sigma)
  maxErr = maxval(abs(zkP - zk) / abs(zk))
  write (*, *) "PBE exchange at zero gradient vs Slater, relative difference ", maxErr
  call check(error, maxErr < 1.0e-13_R64)
  if (allocated(error)) error stop "T_Utils_XcLib: PBE exchange at zero gradient"

  sigma = [1.0e-5_R64, 0.02_R64, 3.0_R64]
  call XcLib_EvaluateGga(zk, vrho, vsigma, funcU, rho, sigma)
  maxErr = 0.0_R64
  do i = 1, 3
    h = 1.0e-4_R64 * sigma(i)
    call XcLib_EvaluateGga(zkP(i:i), vrhoP(i:i), vsigmaP(i:i), funcU, rho(i:i), [sigma(i) + h])
    call XcLib_EvaluateGga(zkM(i:i), vrhoM(i:i), vsigmaM(i:i), funcU, rho(i:i), [sigma(i) - h])
    maxErr = max(maxErr, abs(rho(i) * (zkP(i) - zkM(i)) / (2.0_R64 * h) - vsigma(i)) / abs(vsigma(i)))
  end do
  write (*, *) "PBE exchange v_sigma vs finite difference, relative error ", maxErr
  call check(error, maxErr < 1.0e-6_R64)
  if (allocated(error)) error stop "T_Utils_XcLib: PBE v_sigma"
  call XcLib_End(funcU)

  ! Polarized GGA layout: sigma(3, np), vsigma(3, np)
  call XcLib_Init(funcP, XcLib_FunctionalId("gga_c_pbe"), XCLIB_POLARIZED)
  rho2 = [0.2_R64, 0.1_R64, 0.5_R64, 0.5_R64, 1.0e-3_R64, 2.0e-3_R64]
  sigma3 = [0.01_R64, 0.004_R64, 0.003_R64, 0.1_R64, 0.1_R64, 0.1_R64, 1.0e-6_R64, 1.0e-6_R64, 2.0e-6_R64]
  call XcLib_EvaluateGga(zk, vrho2, vsigma3, funcP, rho2, sigma3)
  call check(error, all(abs(zk) < 1.0_R64) .and. all(zk < 0.0_R64))
  if (allocated(error)) error stop "T_Utils_XcLib: polarized PBE correlation"
  call XcLib_End(funcP)

  !------------------------------------
  ! hybrids
  !------------------------------------

  call XcLib_Init(funcU, XcLib_FunctionalId("hyb_gga_xc_pbeh"), XCLIB_UNPOLARIZED)
  write (*, *) "PBE0 exact exchange fraction ", XcLib_ExactExchangeFraction(funcU)
  call check(error, XcLib_ExactExchangeFraction(funcU), 0.25_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_Utils_XcLib: PBE0 fraction"
  call XcLib_End(funcU)

  call XcLib_Init(funcU, XcLib_FunctionalId("hyb_gga_xc_b3lyp"), XCLIB_UNPOLARIZED)
  call check(error, XcLib_ExactExchangeFraction(funcU), 0.2_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_Utils_XcLib: B3LYP fraction"
  call XcLib_End(funcU)

  call XcLib_Init(funcU, XcLib_FunctionalId("gga_x_pbe"), XCLIB_UNPOLARIZED)
  call check(error, XcLib_ExactExchangeFraction(funcU), 0.0_R64, thr=0.0_R64)
  if (allocated(error)) error stop "T_Utils_XcLib: PBE has no exact exchange"
  call XcLib_End(funcU)

  !------------------------------------
  ! external parameters (1D soft-Coulomb exchange, beta)
  !------------------------------------

  call XcLib_Init(funcU, XcLib_FunctionalId("lda_x_1d_soft"), XCLIB_UNPOLARIZED)
  call check(error, funcU % nExtParams, 1)
  if (allocated(error)) error stop "T_Utils_XcLib: lda_x_1d_soft has one external parameter"
  call check(error, iand(funcU % flags, XCLIB_FLAGS_1D) .ne. 0)
  if (allocated(error)) error stop "T_Utils_XcLib: lda_x_1d_soft is 1D"
  call XcLib_FillExtParamsDefault(params, funcU)
  call check(error, params(1), 1.0_R64, thr=0.0_R64)
  if (allocated(error)) error stop "T_Utils_XcLib: lda_x_1d_soft default beta"
  call XcLib_EvaluateLda(zk, vrho, funcU, rho)
  call XcLib_SetExtParams(funcU, [2.0_R64])
  call XcLib_EvaluateLda(zkP, vrhoP, funcU, rho)
  write (*, *) "lda_x_1d_soft zk(beta = 1), zk(beta = 2) ", zk(2), zkP(2)
  call check(error, abs(zkP(2) - zk(2)) > 1.0e-3_R64)
  if (allocated(error)) error stop "T_Utils_XcLib: lda_x_1d_soft external parameter has no effect"
  call XcLib_End(funcU)

  write (*, *) "T_Utils_XcLib passed"

end program T_Utils_XcLib
