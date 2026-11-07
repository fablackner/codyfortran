!> ctest-requires: libxc
!> SysXc on grid.fem (one center at the origin), density 2|φ|² of the
!> normalized Gaussian φ ∝ e^(-r²/2), i.e. ρ = 2 π^(-3/2) e^(-r²):
!>  - density of a factored 1-RDM (weighted orbital 2φ) and its norm 2;
!>  - Hartree energy 2 sqrt(2/π) (Galerkin Poisson solve of the xc density);
!>  - Slater exchange: energy -(3/4)(3/π)^(1/3) ∫ρ^(4/3) in closed form and
!>    the potential -(3/π)^(1/3) ρ^(1/3) node by node;
!>  - the nodal gradient against the analytic one and the adjoint identity
!>    Σ w u·(G f) = Σ w f (G^† u);
!>  - PBE: v_xc equals the derivative of the discrete E_xc (central finite
!>    difference along a smooth off-center perturbation), unpolarized and
!>    polarized (both spin potentials of an unequal spin density); the
!>    polarized functional of ρ↑ = ρ↓ = ρ/2 equals the unpolarized one;
!>  - Kohn-Sham potential = Hartree + xc; hybrid fraction (PBE0 0.25 and the
!>    exactExchangeFraction override).
program T_SysXc_Fem
  use M_Utils_Types
  use M_Utils_Constants
  use M_Utils_Json
  use M_Grid
  use M_Grid_Fem
  use M_SysInteraction
  use M_SysXc
  use testdrive, only: check, error_type
  implicit none

  type(error_type), allocatable :: error
  real(R64), allocatable :: x(:), y(:), z(:), r2(:)
  complex(R64), allocatable :: phi(:, :), weighted(:, :)
  real(R64), allocatable :: density(:, :), densitySpin(:, :), pert(:), xcPot(:, :), xcPotSpin(:, :), hartreePot(:)
  real(R64), allocatable :: grad(:, :), gradExact(:, :), u(:, :), adj(:), ksPot(:, :)
  real(R64) :: hartreeEnergy, xcEnergy, xcEnergySpin, eExact, maxErr, fd, lin, eps, ePlus, eMinus, lhs, rhs
  real(R64) :: hartreeEnergyKs, xcEnergyKs
  real(R64), allocatable :: scratch2(:, :)
  integer(I32) :: c, s

  call Json_LoadJsonFile(manualFileName_="test/componentTests/SysXc/Fem/T_SysXc_Fem.json", relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysInteraction_Fabricate
  call SysXc_Fabricate
  call Grid_Setup
  call SysInteraction_Setup
  call SysXc_Setup

  write (*, *) "nodes ", Grid_nPoints, ", xc points ", SysXc_nPoints

  call Grid_Fem_FillNodeCoordinate(x, [1.0_R64, 0.0_R64, 0.0_R64])
  call Grid_Fem_FillNodeCoordinate(y, [0.0_R64, 1.0_R64, 0.0_R64])
  call Grid_Fem_FillNodeCoordinate(z, [0.0_R64, 0.0_R64, 1.0_R64])
  r2 = x**2 + y**2 + z**2

  !------------------------------------
  ! density of a factored 1-RDM: ρ = Re conj(φ) (2φ)
  !------------------------------------

  allocate (phi(Grid_nPoints, 1), weighted(Grid_nPoints, 1))
  phi(:, 1) = cmplx(PI**(-0.75_R64) * exp(-0.5_R64 * r2), 0.0_R64, kind=R64) * exp(cmplx(0.0_R64, 0.3_R64 * z, kind=R64))
  weighted = 2.0_R64 * phi

  allocate (density(SysXc_nPoints, 1), xcPot(SysXc_nPoints, 1), hartreePot(SysXc_nPoints))
  density(:, :) = 0.0_R64
  call SysXc_AddDensity(density(:, 1), phi, weighted)
  maxErr = maxval(abs(density(:, 1) - 2.0_R64 * PI**(-1.5_R64) * exp(-r2)))
  write (*, *) "density: max |rho - 2|phi|^2| ", maxErr, ", norm ", SysXc_Integrate(density(:, 1))
  call check(error, maxErr < 1.0e-14_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: density"
  call check(error, SysXc_Integrate(density(:, 1)), 2.0_R64, thr=1.0e-6_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: density norm"

  !------------------------------------
  ! Hartree energy: ½ N² sqrt(2α/π) with N = 2, α = 1
  !------------------------------------

  call SysXc_FillHartreePotential(hartreePot, hartreeEnergy, density(:, 1), 0.0_R64)
  write (*, '(1X, A, F16.12, A, F16.12)') "Hartree energy ", hartreeEnergy, ", exact ", 2.0_R64 * sqrt(2.0_R64 / PI)
  call check(error, hartreeEnergy, 2.0_R64 * sqrt(2.0_R64 / PI), thr=2.0e-6_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: Hartree energy"

  !------------------------------------
  ! Slater exchange in closed form
  !------------------------------------

  call SysXc_FillXcPotential(xcPot, xcEnergy, density)
  ! ∫ρ^(4/3) = N^(4/3) (α/π)² (3π/(4α))^(3/2)
  eExact = -0.75_R64 * (3.0_R64 / PI)**(1.0_R64 / 3.0_R64) * 2.0_R64**(4.0_R64 / 3.0_R64) * PI**(-2) * &
           (0.75_R64 * PI)**1.5_R64
  ! (below Libxc's density threshold the potential is zero)
  maxErr = maxval(abs(xcPot(:, 1) + (3.0_R64 / PI)**(1.0_R64 / 3.0_R64) * density(:, 1)**(1.0_R64 / 3.0_R64)), &
                  mask=density(:, 1) > 1.0e-10_R64)
  write (*, '(1X, A, F16.12, A, F16.12, A, ES10.3)') "Slater exchange energy ", xcEnergy, ", exact ", eExact, &
    ", max potential error ", maxErr
  call check(error, xcEnergy, eExact, thr=1.0e-6_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: Slater exchange energy"
  call check(error, maxErr < 1.0e-12_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: Slater exchange potential"

  !------------------------------------
  ! gradient: accuracy and adjoint identity
  !------------------------------------

  allocate (grad(Grid_nPoints, 3), gradExact(Grid_nPoints, 3), u(Grid_nPoints, 3), adj(Grid_nPoints))
  call SysXc_ApplyGradient(grad, exp(-r2))
  gradExact(:, 1) = -2.0_R64 * x * exp(-r2)
  gradExact(:, 2) = -2.0_R64 * y * exp(-r2)
  gradExact(:, 3) = -2.0_R64 * z * exp(-r2)
  maxErr = maxval(abs(grad - gradExact))
  write (*, *) "gradient of exp(-r^2): max error ", maxErr
  call check(error, maxErr < 5.0e-3_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: gradient accuracy"

  ! Exact for polynomials of degree <= nLocals - 1 per axis away from the
  ! box faces (the Dirichlet values enter the elements touching them)
  call SysXc_ApplyGradient(grad, x**3 * y - 2.0_R64 * z**2 + x * y * z)
  gradExact(:, 1) = 3.0_R64 * x**2 * y + y * z
  gradExact(:, 2) = x**3 + x * z
  gradExact(:, 3) = -4.0_R64 * z + x * y
  maxErr = 0.0_R64
  do c = 1, 3
    maxErr = max(maxErr, maxval(abs(grad(:, c) - gradExact(:, c)), &
                                mask=max(abs(x), abs(y), abs(z)) < 5.5_R64))
  end do
  write (*, *) "gradient of a cubic polynomial (|x|, |y|, |z| < 5.5): max error ", maxErr
  call check(error, maxErr < 1.0e-9_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: gradient is not exact for polynomials"

  u(:, 1) = exp(-0.3_R64 * r2) * (1.0_R64 + y)
  u(:, 2) = exp(-0.5_R64 * ((x - 0.5_R64)**2 + y**2 + z**2))
  u(:, 3) = (1.0_R64 + x) * exp(-0.2_R64 * r2)
  pert = exp(-((x - 0.3_R64)**2 + (y + 0.2_R64)**2 + (z - 0.4_R64)**2))
  call SysXc_ApplyGradient(grad, pert)
  lhs = 0.0_R64
  do c = 1, 3
    lhs = lhs + SysXc_Integrate(u(:, c) * grad(:, c))
  end do
  adj(:) = 0.0_R64
  call SysXc_AddGradientAdjoint(adj, u)
  rhs = SysXc_Integrate(pert * adj)
  write (*, *) "adjoint identity: ", lhs, rhs, ", relative difference ", abs(lhs - rhs) / abs(lhs)
  call check(error, abs(lhs - rhs) / abs(lhs) < 1.0e-12_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: gradient adjoint"

  !------------------------------------
  ! PBE: v_xc = derivative of the discrete E_xc
  !------------------------------------

  call SetFunctionals(["gga_x_pbe", "gga_c_pbe"])
  call check(error, SysXc_gradientQ)
  if (allocated(error)) error stop "T_SysXc_Fem: PBE needs gradients"

  pert = 0.3_R64 * exp(-((x - 0.4_R64)**2 + (y + 0.2_R64)**2 + (z - 0.3_R64)**2) / 0.64_R64) * density(:, 1)
  eps = 1.0e-4_R64

  call SysXc_FillXcPotential(xcPot, xcEnergy, density)
  call SysXc_FillXcPotential(xcPot(:, 1:1), ePlus, reshape(density(:, 1) + eps * pert, [SysXc_nPoints, 1]))
  call SysXc_FillXcPotential(xcPot(:, 1:1), eMinus, reshape(density(:, 1) - eps * pert, [SysXc_nPoints, 1]))
  call SysXc_FillXcPotential(xcPot, xcEnergy, density)
  fd = (ePlus - eMinus) / (2.0_R64 * eps)
  lin = SysXc_Integrate(xcPot(:, 1) * pert)
  write (*, '(1X, A, F16.12, A, 2ES20.12, A, ES10.3)') "PBE: E_xc ", xcEnergy, ", dE finite difference / sum w v drho ", &
    fd, lin, ", relative difference ", abs(fd - lin) / abs(lin)
  call check(error, abs(fd - lin) / abs(lin) < 1.0e-7_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: PBE potential is not the derivative of E_xc"

  ! Polarized functional of equal spin densities = unpolarized
  allocate (densitySpin(SysXc_nPoints, 2), xcPotSpin(SysXc_nPoints, 2))
  densitySpin(:, 1) = 0.5_R64 * density(:, 1)
  densitySpin(:, 2) = 0.5_R64 * density(:, 1)
  call SysXc_FillXcPotential(xcPotSpin, xcEnergySpin, densitySpin)
  maxErr = max(maxval(abs(xcPotSpin(:, 1) - xcPot(:, 1))), maxval(abs(xcPotSpin(:, 2) - xcPot(:, 1))))
  write (*, *) "PBE polarized (equal spins) vs unpolarized: dE ", xcEnergySpin - xcEnergy, ", max dv ", maxErr
  call check(error, abs(xcEnergySpin - xcEnergy) < 1.0e-10_R64 .and. maxErr < 1.0e-8_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: polarized PBE of equal spins"

  ! Unequal spins: both spin potentials are derivatives of E_xc
  densitySpin(:, 1) = 0.65_R64 * density(:, 1)
  densitySpin(:, 2) = 0.35_R64 * density(:, 1) + 0.2_R64 * pert
  call SysXc_FillXcPotential(xcPotSpin, xcEnergySpin, densitySpin)
  allocate (scratch2(SysXc_nPoints, 2))
  do s = 1, 2
    densitySpin(:, s) = densitySpin(:, s) + eps * pert
    call SysXc_FillXcPotential(scratch2, ePlus, densitySpin)
    densitySpin(:, s) = densitySpin(:, s) - 2.0_R64 * eps * pert
    call SysXc_FillXcPotential(scratch2, eMinus, densitySpin)
    densitySpin(:, s) = densitySpin(:, s) + eps * pert
    fd = (ePlus - eMinus) / (2.0_R64 * eps)
    lin = SysXc_Integrate(xcPotSpin(:, s) * pert)
    write (*, '(1X, A, I0, A, 2ES20.12, A, ES10.3)') "PBE polarized, spin ", s, ": finite difference / sum w v drho ", &
      fd, lin, ", relative difference ", abs(fd - lin) / abs(lin)
    call check(error, abs(fd - lin) / abs(lin) < 1.0e-7_R64)
    if (allocated(error)) error stop "T_SysXc_Fem: polarized PBE potential is not the derivative of E_xc"
  end do

  ! Kohn-Sham potential = Hartree + xc
  allocate (ksPot(SysXc_nPoints, 1))
  call SysXc_FillKohnShamPotential(ksPot, hartreeEnergyKs, xcEnergyKs, density, 0.0_R64)
  call SysXc_FillXcPotential(xcPot, xcEnergy, density)
  call SysXc_FillHartreePotential(hartreePot, hartreeEnergy, density(:, 1), 0.0_R64)
  maxErr = maxval(abs(ksPot(:, 1) - hartreePot - xcPot(:, 1)))
  write (*, *) "Kohn-Sham potential - (Hartree + xc): ", maxErr
  call check(error, maxErr < 1.0e-12_R64 .and. abs(hartreeEnergyKs - hartreeEnergy) < 1.0e-13_R64 .and. &
             abs(xcEnergyKs - xcEnergy) < 1.0e-13_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: Kohn-Sham potential"

  !------------------------------------
  ! hybrids
  !------------------------------------

  call SetFunctionals(["hyb_gga_xc_pbeh"])
  write (*, *) "PBE0 exact exchange fraction ", SysXc_exactExchangeFraction
  call check(error, SysXc_exactExchangeFraction, 0.25_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: PBE0 fraction"

  call jsonFile % add("sysXc.fem.libxc.exactExchangeFraction", 0.3_R64)
  call SysXc_Fabricate
  call check(error, SysXc_exactExchangeFraction, 0.3_R64, thr=0.0_R64)
  if (allocated(error)) error stop "T_SysXc_Fem: exactExchangeFraction override"
  call jsonFile % remove("sysXc.fem.libxc.exactExchangeFraction")

  write (*, *) "T_SysXc_Fem passed"

contains

  !> Replaces the functional list and re-fabricates SysXc.
  subroutine SetFunctionals(names)
    character(len=*), intent(in) :: names(:)

    call jsonFile % remove("sysXc.fem.libxc.functionals")
    call jsonFile % add("sysXc.fem.libxc.functionals", names)
    call SysXc_Fabricate
    call SysXc_Setup

  end subroutine

end program T_SysXc_Fem
