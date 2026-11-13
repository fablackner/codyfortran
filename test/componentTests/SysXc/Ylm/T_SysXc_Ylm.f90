!> ctest-requires: libxc
!> SysXc on a Ylm grid (lmax 3, mmax 2: general mode with nPhi = 9) for two
!> complex orbitals of mixed l and m (non-spherical density):
!>  - density integral = Σ ⟨φ|φ⟩ and the factored form Re conj(φ) (2φ);
!>  - Hartree energy and J·φ from the xc density equal the Gaunt-product path
!>    of SysInteraction (SysInteraction_FillDirectPotential,
!>    MultiplyWithInteractionPotential) to round-off (the angular quadrature
!>    is exact for these band-limited fields);
!>  - the local potential acts Hermitian in the Ylm metric;
!>  - gradient of e^(-r²)(1 + x + z/2) in the spherical frame against the
!>    analytic components, adjoint identity Σ w u·(G f) = Σ w f (G^† u);
!>  - PBE: v_xc = derivative of the discrete E_xc (unpolarized, polarized);
!>  - axial mode (nPhi = 1) and general mode agree for an axial density.
program T_SysXc_Ylm
  use M_Utils_Types
  use M_Utils_Constants
  use M_Utils_Json
  use M_Grid
  use M_Grid_Ylm
  use M_SysInteraction
  use M_SysXc
  use M_SysXc_Ylm
  use testdrive, only: check, error_type
  implicit none

  type(error_type), allocatable :: error
  complex(R64), allocatable :: orbs(:, :), orbsAxial(:, :), orbTmp(:), orbTmp2(:)
  complex(R64), allocatable :: potGeneric(:)
  real(R64), allocatable :: r(:), density(:, :), density2(:), hartreePot(:), xcPot(:, :), densitySpin(:, :)
  real(R64), allocatable :: pert(:), f(:), grad(:, :), gradExact(:, :), u(:, :), adj(:), scratch(:, :)
  real(R64), allocatable :: xP(:), zP(:), rP(:), cosP(:), sinP(:), phiP(:)
  real(R64) :: norm, hartreeEnergy, hartreeGeneric, xcEnergy, xcGeneral, xcAxial, maxErr, lhs, rhs
  real(R64) :: eps, ePlus, eMinus, fd, lin
  complex(R64) :: avb, bva
  integer(I32) :: nRad, nT, nPh, i, q, k, p, c, s

  call Json_LoadJsonFile(manualFileName_="test/componentTests/SysXc/Ylm/T_SysXc_Ylm.json", relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysInteraction_Fabricate
  call SysXc_Fabricate
  call Grid_Setup
  call SysInteraction_Setup
  call SysXc_Setup

  nRad = Grid_Ylm_nRadial
  nT = SysXc_Ylm_nTheta
  nPh = SysXc_Ylm_nPhi
  r = Grid_Ylm_radialPoints
  write (*, *) "radial points ", nRad, ", nTheta ", nT, ", nPhi ", nPh, ", xc points ", SysXc_nPoints
  call check(error,.not. SysXc_Ylm_axialQ .and. nPh .eq. 9)
  if (allocated(error)) error stop "T_SysXc_Ylm: general mode with nPhi = 4 mmax + 1 expected"

  !------------------------------------
  ! orbitals of mixed l and m
  !------------------------------------

  allocate (orbs(Grid_nPoints, 2), orbTmp(Grid_nPoints), orbTmp2(Grid_nPoints))
  orbs(:, :) = (0.0_R64, 0.0_R64)
  call Grid_Ylm_SetLmComponent(orbs(:, 1), 0, 0, cmplx(exp(-0.5_R64 * r**2), 0.0_R64, kind=R64))
  call Grid_Ylm_SetLmComponent(orbs(:, 1), 1, 0, cmplx(0.3_R64 * r * exp(-0.5_R64 * r**2), 0.0_R64, kind=R64))
  call Grid_Ylm_SetLmComponent(orbs(:, 1), 2, 1, (0.1_R64, 0.05_R64) * r**2 * exp(-0.5_R64 * r**2))
  call Grid_Ylm_SetLmComponent(orbs(:, 2), 1, 1, cmplx(r * exp(-0.8_R64 * r**2), 0.0_R64, kind=R64))
  call Grid_Ylm_SetLmComponent(orbs(:, 2), 1, -1, (0.0_R64, 0.5_R64) * r * exp(-0.8_R64 * r**2))
  call Grid_Ylm_SetLmComponent(orbs(:, 2), 0, 0, cmplx(0.2_R64 * exp(-r**2), 0.0_R64, kind=R64))
  call Grid_Ylm_SetLmComponent(orbs(:, 2), 3, -2, cmplx(0.05_R64 * r**3 * exp(-0.8_R64 * r**2), 0.0_R64, kind=R64))

  !------------------------------------
  ! density
  !------------------------------------

  allocate (density(SysXc_nPoints, 1), density2(SysXc_nPoints))
  density(:, :) = 0.0_R64
  call SysXc_AddDensity(density(:, 1), orbs, orbs)
  norm = real(Grid_InnerProduct(orbs(:, 1), orbs(:, 1)) + Grid_InnerProduct(orbs(:, 2), orbs(:, 2)), kind=R64)
  write (*, *) "density integral ", SysXc_Integrate(density(:, 1)), ", sum of norms ", norm
  call check(error, abs(SysXc_Integrate(density(:, 1)) - norm) < 1.0e-13_R64 * norm)
  if (allocated(error)) error stop "T_SysXc_Ylm: density integral"
  call check(error, minval(density(:, 1)) >= 0.0_R64)
  if (allocated(error)) error stop "T_SysXc_Ylm: negative density"

  density2(:) = 0.0_R64
  call SysXc_AddDensity(density2, orbs, 2.0_R64 * orbs)
  maxErr = maxval(abs(density2 - 2.0_R64 * density(:, 1)))
  call check(error, maxErr < 1.0e-14_R64 * maxval(density))
  if (allocated(error)) error stop "T_SysXc_Ylm: factored density"

  !------------------------------------
  ! Hartree: xc density vs the Gaunt-product path
  !------------------------------------

  allocate (hartreePot(SysXc_nPoints))
  call SysXc_FillHartreePotential(hartreePot, hartreeEnergy, density(:, 1), 0.0_R64)
  call SysInteraction_FillDirectPotential(potGeneric, orbs, orbs, 0.0_R64)
  hartreeGeneric = 0.0_R64
  do i = 1, 2
    call SysInteraction_MultiplyWithInteractionPotential(orbTmp, potGeneric, orbs(:, i))
    hartreeGeneric = hartreeGeneric + 0.5_R64 * real(Grid_InnerProduct(orbs(:, i), orbTmp), kind=R64)
  end do
  write (*, '(1X, A, F18.14, A, F18.14)') "Hartree energy: xc points ", hartreeEnergy, ", Gaunt path ", hartreeGeneric
  call check(error, abs(hartreeEnergy - hartreeGeneric) < 1.0e-12_R64 * abs(hartreeGeneric))
  if (allocated(error)) error stop "T_SysXc_Ylm: Hartree energy"

  maxErr = 0.0_R64
  do i = 1, 2
    call SysXc_MultiplyWithPotential(orbTmp, hartreePot, orbs(:, i))
    call SysInteraction_MultiplyWithInteractionPotential(orbTmp2, potGeneric, orbs(:, i))
    maxErr = max(maxErr, maxval(abs(orbTmp - orbTmp2)) / maxval(abs(orbTmp2)))
  end do
  write (*, *) "J phi: xc points vs Gaunt path, max relative difference ", maxErr
  call check(error, maxErr < 1.0e-11_R64)
  if (allocated(error)) error stop "T_SysXc_Ylm: Hartree potential application"

  !------------------------------------
  ! LDA potential: Hermitian
  !------------------------------------

  allocate (xcPot(SysXc_nPoints, 1))
  call SysXc_FillXcPotential(xcPot, xcEnergy, density)
  call SysXc_MultiplyWithPotential(orbTmp, xcPot(:, 1), orbs(:, 2))
  avb = Grid_InnerProduct(orbs(:, 1), orbTmp)
  call SysXc_MultiplyWithPotential(orbTmp, xcPot(:, 1), orbs(:, 1))
  bva = Grid_InnerProduct(orbs(:, 2), orbTmp)
  write (*, *) "LDA: E_xc ", xcEnergy, ", <a|v b> - conj <b|v a> ", abs(avb - conjg(bva))
  call check(error, abs(avb - conjg(bva)) < 1.0e-13_R64 * abs(avb))
  if (allocated(error)) error stop "T_SysXc_Ylm: local potential not Hermitian"

  !------------------------------------
  ! gradient in the spherical frame
  !------------------------------------

  allocate (rP(SysXc_nPoints), cosP(SysXc_nPoints), sinP(SysXc_nPoints), phiP(SysXc_nPoints))
  do k = 1, nPh
    do q = 1, nT
      do i = 1, nRad
        p = i + nRad * (q - 1) + nRad * nT * (k - 1)
        rP(p) = r(i)
        cosP(p) = SysXc_Ylm_cosTheta(q)
        sinP(p) = sqrt(1.0_R64 - SysXc_Ylm_cosTheta(q)**2)
        phiP(p) = 2.0_R64 * PI * (k - 1) / real(nPh, kind=R64)
      end do
    end do
  end do
  xP = rP * sinP * cos(phiP)
  zP = rP * cosP

  allocate (grad(SysXc_nPoints, 3), gradExact(SysXc_nPoints, 3), u(SysXc_nPoints, 3), adj(SysXc_nPoints))
  f = exp(-rP**2) * (1.0_R64 + xP + 0.5_R64 * zP)
  call SysXc_ApplyGradient(grad, f)
  gradExact(:, 1) = exp(-rP**2) * (-2.0_R64 * rP * (1.0_R64 + xP + 0.5_R64 * zP) + sinP * cos(phiP) + 0.5_R64 * cosP)
  gradExact(:, 2) = exp(-rP**2) * (cosP * cos(phiP) - 0.5_R64 * sinP)
  gradExact(:, 3) = -exp(-rP**2) * sin(phiP)
  do c = 1, 3
    write (*, '(1X, A, I0, A, ES10.3)') "gradient component ", c, ": max error ", maxval(abs(grad(:, c) - gradExact(:, c)))
  end do
  call check(error, maxval(abs(grad(:, 2:3) - gradExact(:, 2:3))) < 1.0e-12_R64)
  if (allocated(error)) error stop "T_SysXc_Ylm: angular gradient"
  p = maxloc(abs(grad(:, 1) - gradExact(:, 1)), dim=1)
  write (*, '(1X, A, F10.5)') "  (largest radial error at r = ", rP(p)
  call check(error, maxval(abs(grad(:, 1) - gradExact(:, 1))) < 1.0e-3_R64)
  if (allocated(error)) error stop "T_SysXc_Ylm: radial gradient"

  ! Exact for r f polynomial of degree <= nLocals - 1 (f = 1 + r/2 - r^3/10)
  f = (1.0_R64 + 0.5_R64 * rP - 0.1_R64 * rP**3) * (1.0_R64 + xP)
  call SysXc_ApplyGradient(grad, f)
  gradExact(:, 1) = (0.5_R64 - 0.3_R64 * rP**2) * (1.0_R64 + xP) + (1.0_R64 + 0.5_R64 * rP - 0.1_R64 * rP**3) * sinP * cos(phiP)
  maxErr = maxval(abs(grad(:, 1) - gradExact(:, 1))) / maxval(abs(gradExact(:, 1)))
  write (*, '(1X, A, ES10.3)') "radial gradient of a polynomial: max relative error ", maxErr
  call check(error, maxErr < 1.0e-11_R64)
  if (allocated(error)) error stop "T_SysXc_Ylm: radial gradient is not exact for polynomials"

  u(:, 1) = exp(-0.3_R64 * rP**2) * (1.0_R64 + zP)
  u(:, 2) = exp(-0.5_R64 * rP**2) * (0.5_R64 + xP)
  u(:, 3) = exp(-0.4_R64 * rP**2) * (1.0_R64 - 0.3_R64 * xP * zP)
  call SysXc_ApplyGradient(grad, density(:, 1))
  lhs = 0.0_R64
  do c = 1, 3
    lhs = lhs + SysXc_Integrate(u(:, c) * grad(:, c))
  end do
  adj(:) = 0.0_R64
  call SysXc_AddGradientAdjoint(adj, u)
  rhs = SysXc_Integrate(density(:, 1) * adj)
  write (*, *) "adjoint identity: ", lhs, rhs, ", relative difference ", abs(lhs - rhs) / abs(lhs)
  call check(error, abs(lhs - rhs) < 1.0e-12_R64 * abs(lhs))
  if (allocated(error)) error stop "T_SysXc_Ylm: gradient adjoint"

  !------------------------------------
  ! PBE: v_xc = derivative of the discrete E_xc
  !------------------------------------

  call SetFunctionals(["gga_x_pbe", "gga_c_pbe"])
  pert = 0.3_R64 * exp(-((xP - 0.4_R64)**2 + (zP + 0.3_R64)**2 + (rP * sinP * sin(phiP))**2)) * density(:, 1)
  eps = 1.0e-4_R64
  allocate (scratch(SysXc_nPoints, 1))
  call SysXc_FillXcPotential(scratch, ePlus, reshape(density(:, 1) + eps * pert, [SysXc_nPoints, 1]))
  call SysXc_FillXcPotential(scratch, eMinus, reshape(density(:, 1) - eps * pert, [SysXc_nPoints, 1]))
  call SysXc_FillXcPotential(xcPot, xcEnergy, density)
  fd = (ePlus - eMinus) / (2.0_R64 * eps)
  lin = SysXc_Integrate(xcPot(:, 1) * pert)
  write (*, '(1X, A, F16.12, A, 2ES20.12, A, ES10.3)') "PBE: E_xc ", xcEnergy, ", finite difference / sum w v drho ", &
    fd, lin, ", relative difference ", abs(fd - lin) / abs(lin)
  call check(error, abs(fd - lin) / abs(lin) < 1.0e-7_R64)
  if (allocated(error)) error stop "T_SysXc_Ylm: PBE potential is not the derivative of E_xc"

  allocate (densitySpin(SysXc_nPoints, 2))
  deallocate (scratch)
  allocate (scratch(SysXc_nPoints, 2))
  densitySpin(:, 1) = 0.7_R64 * density(:, 1)
  densitySpin(:, 2) = 0.3_R64 * density(:, 1) + 0.1_R64 * pert
  call SysXc_FillXcPotential(scratch, xcEnergy, densitySpin)
  do s = 1, 2
    xcPot(:, 1) = scratch(:, s)
    densitySpin(:, s) = densitySpin(:, s) + eps * pert
    call SysXc_FillXcPotential(scratch, ePlus, densitySpin)
    densitySpin(:, s) = densitySpin(:, s) - 2.0_R64 * eps * pert
    call SysXc_FillXcPotential(scratch, eMinus, densitySpin)
    densitySpin(:, s) = densitySpin(:, s) + eps * pert
    call SysXc_FillXcPotential(scratch, xcEnergy, densitySpin)
    fd = (ePlus - eMinus) / (2.0_R64 * eps)
    lin = SysXc_Integrate(xcPot(:, 1) * pert)
    write (*, '(1X, A, I0, A, 2ES20.12, A, ES10.3)') "PBE polarized, spin ", s, ": finite difference / sum w v drho ", &
      fd, lin, ", relative difference ", abs(fd - lin) / abs(lin)
    call check(error, abs(fd - lin) / abs(lin) < 1.0e-7_R64)
    if (allocated(error)) error stop "T_SysXc_Ylm: polarized PBE potential is not the derivative of E_xc"
  end do

  !------------------------------------
  ! axial mode = general mode for an axial density (orbitals of definite m)
  !------------------------------------

  allocate (orbsAxial(Grid_nPoints, 2))
  orbsAxial(:, :) = (0.0_R64, 0.0_R64)
  call Grid_Ylm_SetLmComponent(orbsAxial(:, 1), 0, 0, cmplx(exp(-0.5_R64 * r**2), 0.0_R64, kind=R64))
  call Grid_Ylm_SetLmComponent(orbsAxial(:, 1), 1, 0, cmplx(0.3_R64 * r * exp(-0.5_R64 * r**2), 0.0_R64, kind=R64))
  call Grid_Ylm_SetLmComponent(orbsAxial(:, 2), 1, 1, cmplx(r * exp(-0.8_R64 * r**2), 0.0_R64, kind=R64))
  call Grid_Ylm_SetLmComponent(orbsAxial(:, 2), 2, 1, (0.1_R64, 0.2_R64) * r**2 * exp(-0.8_R64 * r**2))

  xcGeneral = AxialTestEnergy(orbsAxial)
  call jsonFile % add("sysXc.ylm.nPhi", 1)
  call SysXc_Fabricate
  call SysXc_Setup
  call check(error, SysXc_Ylm_axialQ)
  if (allocated(error)) error stop "T_SysXc_Ylm: nPhi = 1 must select the axial mode"
  xcAxial = AxialTestEnergy(orbsAxial)
  write (*, '(1X, A, F18.14, A, F18.14)') "PBE E_xc of an axial density: general ", xcGeneral, ", axial ", xcAxial
  call check(error, abs(xcAxial - xcGeneral) < 1.0e-12_R64 * abs(xcGeneral))
  if (allocated(error)) error stop "T_SysXc_Ylm: axial and general mode differ"

  write (*, *) "T_SysXc_Ylm passed"

contains

  !> Replaces the functional list and re-fabricates SysXc.
  subroutine SetFunctionals(names)
    character(len=*), intent(in) :: names(:)

    call jsonFile % remove("sysXc.ylm.libxc.functionals")
    call jsonFile % add("sysXc.ylm.libxc.functionals", names)
    call SysXc_Fabricate
    call SysXc_Setup

  end subroutine

  !> E_xc (current functional and mode) of the density of the given orbitals.
  function AxialTestEnergy(orbitals) result(res)
    complex(R64), intent(in), contiguous :: orbitals(:, :)
    real(R64) :: res

    real(R64), allocatable :: rho(:, :), v(:, :)

    allocate (rho(SysXc_nPoints, 1), v(SysXc_nPoints, 1))
    rho(:, :) = 0.0_R64
    call SysXc_AddDensity(rho(:, 1), orbitals, orbitals)
    call SysXc_FillXcPotential(v, res, rho)

  end function

end program T_SysXc_Ylm
