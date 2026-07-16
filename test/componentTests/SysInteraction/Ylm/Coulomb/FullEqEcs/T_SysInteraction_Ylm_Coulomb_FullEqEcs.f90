program T_SysInteraction_Ylm_Coulomb_FullEqEcs
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Ylm
  use M_Grid_Ylm_FedvrEcs
  use M_SysInteraction
  use M_SysInteraction_Ylm
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: i, l
  real(R64) :: avgError, avgPot, lambda
  complex(R64) :: z, analyticalPot
  complex(R64), allocatable :: srcLm(:), potLm(:)
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/SysInteraction/Ylm/Coulomb/FullEqEcs/T_SysInteraction_Ylm_Coulomb_FullEqEcs.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysInteraction_Fabricate

  call Grid_Setup
  call SysInteraction_Setup

  ! Allocate arrays
  allocate (srcLm(Grid_Ylm_nRadial), potLm(Grid_Ylm_nRadial))

  !==========================================
  call Say_Section("perform action")
  !==========================================

  ! --- l = 0: normalized exponential density ρ(z) = (λ³/8π) e^{-λz} whose
  ! potential V(z) = [1 - e^{-λz}(1 + λz/2)]/z carries the physical Coulomb
  ! tail 1/z, which is genuinely complex in the ECS region
  l = 0
  lambda = 3.0_R64
  do i = 1, Grid_Ylm_nRadial
    z = Grid_Ylm_FedvrEcs_contourPoints(i)
    srcLm(i) = (lambda**3 / (8.0_R64 * PI)) * exp(-lambda * z)
    srcLm(i) = srcLm(i) * Grid_Ylm_radialMetricWeights(i)
  end do

  call SysInteraction_Ylm_FillInteractionPotentialRadial(potLm, srcLm, l, 0, 0.0_R64)

  avgError = 0.0_R64
  avgPot = 0.0_R64
  do i = 1, Grid_Ylm_nRadial
    z = Grid_Ylm_FedvrEcs_contourPoints(i)
    analyticalPot = (1.0_R64 - exp(-lambda * z) * (1.0_R64 + 0.5_R64 * lambda * z)) / z
    avgError = avgError + abs(potLm(i) - analyticalPot)
    avgPot = avgPot + abs(analyticalPot)
  end do
  write (*, '(A, E12.5)') "Average error (l=0, Coulomb tail): ", avgError / avgPot

  call check(error, avgError / avgPot < 1e-7_R64)
  if (allocated(error)) error stop "T_SysInteraction_Ylm_Coulomb_FullEqEcs l=0 failure"

  ! --- l = 1: exponentially confined manufactured solution V(z) = z e^{-λz}
  ! exercising the centrifugal term (same pair as the FullEq test, continued
  ! onto the contour)
  l = 1
  lambda = 5.0_R64
  do i = 1, Grid_Ylm_nRadial
    z = Grid_Ylm_FedvrEcs_contourPoints(i)
    srcLm(i) = -(lambda / (4.0_R64 * PI)) * (lambda * z - 2.0_R64 * (l + 1)) * z**(l - 1) * exp(-lambda * z)
    srcLm(i) = srcLm(i) * Grid_Ylm_radialMetricWeights(i)
  end do

  call SysInteraction_Ylm_FillInteractionPotentialRadial(potLm, srcLm, l, 0, 0.0_R64)

  avgError = 0.0_R64
  avgPot = 0.0_R64
  do i = 1, Grid_Ylm_nRadial
    z = Grid_Ylm_FedvrEcs_contourPoints(i)
    analyticalPot = z**l * exp(-lambda * z)
    avgError = avgError + abs(potLm(i) - analyticalPot)
    avgPot = avgPot + abs(analyticalPot)
  end do
  write (*, '(A, E12.5)') "Average error (l=1, centrifugal): ", avgError / avgPot

  call check(error, avgError / avgPot < 1e-7_R64)
  if (allocated(error)) error stop "T_SysInteraction_Ylm_Coulomb_FullEqEcs l=1 failure"

  ! Deallocate all arrays
  deallocate (srcLm, potLm)

  !==========================================
  call Say_Goodbye
  !==========================================

end program
