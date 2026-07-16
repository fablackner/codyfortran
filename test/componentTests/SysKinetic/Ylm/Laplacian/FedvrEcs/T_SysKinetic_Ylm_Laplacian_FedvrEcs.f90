program T_SysKinetic_Ylm_Laplacian_FedvrEcs
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Ylm
  use M_Grid_Ylm_FedvrEcs
  use M_SysKinetic
  use M_SysKinetic_Ylm
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: i, l
  real(R64) :: avgError, avgResult, lambda, a
  complex(R64) :: z, zEnd, analyticalResult, cnorm, cnormExact
  complex(R64), allocatable :: wavefunc(:), kineticResult(:), fFull(:)
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/SysKinetic/Ylm/Laplacian/FedvrEcs/T_SysKinetic_Ylm_Laplacian_FedvrEcs.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate

  call Grid_Setup
  call SysKinetic_Setup

  ! Allocate only needed arrays
  allocate (wavefunc(Grid_Ylm_nRadial), kineticResult(Grid_Ylm_nRadial))

  l = 1
  lambda = 1.0_R64

  ! Initialize test wavefunction f(z)=z^l e^{-λ z} on the ECS contour
  do i = 1, Grid_Ylm_nRadial
    z = Grid_Ylm_FedvrEcs_contourPoints(i)
    wavefunc(i) = z**l * exp(-lambda * z)
  end do

  ! Apply kinetic operator T = -1/2 ∇² with ∇² f = (1/z) d²(z f)/dz² - l(l+1)/z² f
  call SysKinetic_Ylm_MultiplyWithRadialKineticOp(kineticResult, wavefunc, l, 0, 0.0_R64)

  ! Analytical kinetic result (analytic continuation onto the contour):
  ! T f = -0.5 * λ * z^{l-1} e^{-λ z} (λ z - 2(l+1))
  ! The included free endpoint at rmax carries the weak-form boundary term
  ! (natural boundary condition), so it is excluded from the error metric
  avgError = 0.0_R64
  avgResult = 0.0_R64
  do i = 1, Grid_Ylm_nRadial - 1
    z = Grid_Ylm_FedvrEcs_contourPoints(i)
    analyticalResult = -0.5_R64 * lambda * z**(l - 1) * exp(-lambda * z) * (lambda * z - 2.0_R64 * (l + 1))
    avgError = avgError + abs(kineticResult(i) - analyticalResult)
    avgResult = avgResult + abs(analyticalResult)
  end do
  write (*, '(A, 3E12.5)') "Average error (Kinetic): ", avgError / avgResult

  ! Check if the average error is within the tolerance
  call check(error, avgError / avgResult < 1e-5_R64)
  if (allocated(error)) error stop "SysKinetic_Ylm_Laplacian_FedvrEcs kinetic failure"

  ! Verify the c-product inner product bound by the ECS grid:
  ! c-norm of f in a single (l,m) channel is the contour integral of f² z²,
  ! for f = z e^{-λ z} the antiderivative of z⁴ e^{-a z} with a = 2λ gives
  ! ∫₀^{zEnd} = 24/a⁵ - e^{-a zEnd}(zEnd⁴/a + 4 zEnd³/a² + 12 zEnd²/a³ + 24 zEnd/a⁴ + 24/a⁵)
  allocate (fFull(Grid_nPoints))
  fFull = (0.0_R64, 0.0_R64)
  call Grid_Ylm_AddLmComponent(fFull, l, 0, wavefunc)

  cnorm = Grid_InnerProduct(fFull, fFull)

  a = 2.0_R64 * lambda
  zEnd = Grid_Ylm_FedvrEcs_contourPoints(Grid_Ylm_nRadial)
  cnormExact = 24.0_R64 / a**5 - exp(-a * zEnd) * &
               (zEnd**4 / a + 4.0_R64 * zEnd**3 / a**2 + 12.0_R64 * zEnd**2 / a**3 + &
                24.0_R64 * zEnd / a**4 + 24.0_R64 / a**5)
  write (*, '(A, 3E12.5)') "Relative error (c-norm): ", abs(cnorm - cnormExact) / abs(cnormExact)

  ! Check if the c-norm matches the analytic contour integral
  call check(error, abs(cnorm - cnormExact) / abs(cnormExact) < 1e-6_R64)
  if (allocated(error)) error stop "SysKinetic_Ylm_Laplacian_FedvrEcs c-norm failure"

  deallocate (wavefunc, kineticResult, fFull)

  !==========================================
  call Say_Goodbye
  !==========================================

end program
