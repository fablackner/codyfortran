program T_SysInteraction_Ylm_Coulomb_FullEq
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Ylm
  use M_SysInteraction
  use M_SysInteraction_Ylm
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: i, l
  real(R64) :: r, analyticalPot, avgError, avgPot
  complex(R64), allocatable :: srcLm(:), potLm(:)
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  real(R64) :: lambda

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/T_SysInteraction_Ylm_Coulomb_FullEq.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysInteraction_Fabricate

  call Grid_Setup
  call SysInteraction_Setup

  ! Allocate arrays
  allocate (srcLm(Grid_Ylm_nRadial), potLm(Grid_Ylm_nRadial))

  ! Use a single l value
  l = 1

  ! Initialize with a std charge density with a known analytical solution
  lambda = 5.0_R64  ! or any positive number
  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    srcLm(i) = -(lambda / (4.0_R64 * PI)) * (lambda * r - 2.0_R64 * (l + 1)) * r**(l - 1) * exp(-lambda * r)
    srcLm(i) = srcLm(i) * Grid_Ylm_radialWeights(i)
  end do

  !==========================================
  call Say_Section("perform action")
  !==========================================

  ! Solve using PoissonEqStd
  call SysInteraction_Ylm_FillInteractionPotentialRadial(potLm, srcLm, l, 0, 0.0_R64)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  ! Save potential to file and calculate average error
  avgError = 0.0_R64
  avgPot = 0.0_R64
  do i = 1, Grid_Ylm_nRadial
    r = Grid_Ylm_radialPoints(i)
    analyticalPot = r**l * exp(-lambda * r)
    avgError = avgError + abs(real(potLm(i)) - analyticalPot)
    avgPot = avgPot + analyticalPot
  end do
  write (*, '(A, E12.5)') "Average error: ", avgError / avgPot

  ! Check if the average error is within the tolerance
  call check(error, avgError / avgPot < 1e-7_R64)
  if (allocated(error)) error stop "T_SysInteraction_02_Ylm_Coulomb_FullEq failure"

  ! Deallocate all arrays
  deallocate (srcLm, potLm)

  !==========================================
  call Say_Goodbye
  !==========================================

end program
