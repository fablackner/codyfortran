!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Radial Poisson solve of a manufactured pair, shared by the tests of the Ylm
! Coulomb backends on real FEDVR grids (stdImpl, twoScan, fullEq, blockEq);
! included after `contains`. The ECS backend has its own test (complex
! contour, l = 0 Coulomb tail).
!
! For l = 1 and lambda = 5 the source
!   rho(r) = -(lambda / 4 pi) (lambda r - 2 (l + 1)) r^(l - 1) exp(-lambda r)
! has the potential V(r) = r^l exp(-lambda r), which has decayed long before
! rmax = 10, so the boundary condition of a solver does not enter.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Loads jsonFileName, solves the radial Poisson problem of the manufactured
!> pair with SysInteraction_Ylm_FillInteractionPotentialRadial and stops if
!> the average relative error of the potential,
!> sum_i |Re V_num(r_i) - V(r_i)| / sum_i |V(r_i)|, reaches tolerance.
subroutine CheckManufacturedPoisson(jsonFileName, tolerance)
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Ylm
  use M_SysInteraction
  use M_SysInteraction_Ylm
  use testdrive, only: check, error_type

  character(len=*), intent(in) :: jsonFileName
  real(R64), intent(in) :: tolerance

  integer(I32), parameter :: l = 1
  real(R64), parameter :: lambda = 5.0_R64

  integer(I32) :: i
  real(R64) :: r, analyticalPot, avgError, avgPot
  complex(R64), allocatable :: srcLm(:), potLm(:)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysInteraction_Fabricate

  call Grid_Setup
  call SysInteraction_Setup

  ! The source density times the radial quadrature weights
  allocate (srcLm(Grid_Ylm_nRadial), potLm(Grid_Ylm_nRadial))
  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    srcLm(i) = -(lambda / (4.0_R64 * PI)) * (lambda * r - 2.0_R64 * (l + 1)) * r**(l - 1) * exp(-lambda * r)
    srcLm(i) = srcLm(i) * Grid_Ylm_radialWeights(i)
  end do

  !==========================================
  call Say_Section("perform action")
  !==========================================

  call SysInteraction_Ylm_FillInteractionPotentialRadial(potLm, srcLm, l, 0, 0.0_R64)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  avgError = 0.0_R64
  avgPot = 0.0_R64
  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    analyticalPot = r**l * exp(-lambda * r)
    avgError = avgError + abs(real(potLm(i)) - analyticalPot)
    avgPot = avgPot + analyticalPot
  end do
  write (*, '(A, E12.5)') "Average error: ", avgError / avgPot

  call check(error, avgError / avgPot < tolerance)
  if (allocated(error)) error stop "CheckManufacturedPoisson: average error of the potential above the tolerance"

end subroutine
