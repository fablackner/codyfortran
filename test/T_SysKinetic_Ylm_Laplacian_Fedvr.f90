program T_SysKinetic_Ylm_Laplacian_Fedvr
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Ylm
  use M_SysInteraction
  use M_SysInteraction_Ylm
  use M_SysKinetic
  use M_SysKinetic_Ylm
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: i, l
  real(R64) :: r, avgError, avgResult
  complex(R64), allocatable :: wavefunc(:), kineticResult(:)
  real(R64) :: analyticalResult
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  real(R64) :: lambda

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/T_SysKinetic_Ylm_Laplacian_Fedvr.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate

  call Grid_Setup
  call SysKinetic_Setup

  ! Allocate only needed arrays
  allocate (wavefunc(Grid_Ylm_nRadial), kineticResult(Grid_Ylm_nRadial))

  l = 1
  lambda = 5.0_R64

  ! Initialize test wavefunction f(r)=r^l e^{-λ r}
  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    wavefunc(i) = r**l * exp(-lambda * r)
  end do

  ! Apply kinetic operator T = -1/2 ∇² with ∇² f = (1/r) d²(r f)/dr² - l(l+1)/r² f
  call SysKinetic_Ylm_MultiplyWithRadialKineticOp(kineticResult, wavefunc, l, 0, 0.0_R64)

  ! Analytical kinetic result:
  ! T f = -0.5 * [ (1/r) d²(r f)/dr² - l(l+1)/r² f ]
  ! For f = r^l e^{-λ r} this simplifies to:
  ! T f = -0.5 * λ * r^{l-1} e^{-λ r} (λ r - 2(l+1))
  avgError = 0.0_R64
  avgResult = 0.0_R64
  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    analyticalResult = -0.5_R64 * lambda * r**(l - 1) * exp(-lambda * r) * (lambda * r - 2.0_R64 * (l + 1))
    avgError = avgError + abs(real(kineticResult(i)) - analyticalResult)
    avgResult = avgResult + analyticalResult
  end do
  write (*, '(A, 3E12.5)') "Average error (Kinetic): ", avgError / avgResult

  ! Check if the average error is within the tolerance
  call check(error, avgError / avgResult < 1e-5_R64)
  if (allocated(error)) error stop "SysKinetic_Ylm_Laplacian_Fedvr failure"

  deallocate (wavefunc, kineticResult)

  !==========================================
  call Say_Goodbye
  !==========================================

end program
