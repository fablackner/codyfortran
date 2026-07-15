program T_SysInteraction_Prolate_Coulomb_StdImpl
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Prolate
  use M_SysInteraction
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: iGrid
  real(R64) :: sigma, r, analyticalPot, avgError, avgPot, hartreeAnalytic
  complex(R64) :: hartree
  complex(R64), allocatable :: orb(:), src(:), pot(:), dOrb(:)
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/SysInteraction/Prolate/Coulomb/StdImpl/T_SysInteraction_Prolate_Coulomb_StdImpl.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysInteraction_Fabricate

  call Grid_Setup
  call SysInteraction_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  ! Orbital whose density is a spherical Gaussian centered at the bond
  ! midpoint: |psi|^2 = exp(-r^2/sigma^2) with r^2 = a^2 (xi^2 + eta^2 - 1).
  ! Its Coulomb potential is analytic:
  !   V(r) = pi^(3/2) sigma^3 erf(r/sigma) / r
  ! and the Hartree integral is
  !   <psi | V psi> = sqrt(2) pi^(5/2) sigma^5.
  sigma = 1.0_R64

  allocate (orb(Grid_nPoints), dOrb(Grid_nPoints))
  do iGrid = 1, Grid_nPoints
    r = Grid_Prolate_a * sqrt(Grid_Prolate_xiCoord(iGrid)**2 + Grid_Prolate_etaCoord(iGrid)**2 - 1.0_R64)
    orb(iGrid) = sqrt(TWOPI) * exp(-r**2 / (2.0_R64 * sigma**2))
  end do

  call SysInteraction_FillInteractionSrc(src, orb, orb)
  call SysInteraction_FillInteractionPotential(pot, src, 0.0_R64)

  !------------------------------------
  ! pointwise potential vs. analytic erf solution (m = 0 channel stores
  ! sqrt(2 pi) V)
  !------------------------------------

  avgError = 0.0_R64
  avgPot = 0.0_R64
  do iGrid = 1, Grid_Prolate_nSpatial
    r = Grid_Prolate_a * sqrt(Grid_Prolate_xiCoord(iGrid)**2 + Grid_Prolate_etaCoord(iGrid)**2 - 1.0_R64)
    analyticalPot = PI**1.5_R64 * sigma**3 * erf(r / sigma) / r
    avgError = avgError + abs(real(pot(iGrid), kind=R64) / sqrt(TWOPI) - analyticalPot)
    avgPot = avgPot + analyticalPot
  end do
  write (*, '(A, E12.5)') "Average potential error: ", avgError / avgPot

  call check(error, avgError / avgPot < 1e-8_R64)
  if (allocated(error)) error stop "T_SysInteraction_Prolate potential failure"

  !------------------------------------
  ! Hartree integral via MultiplyWithInteractionPotential
  !------------------------------------

  call SysInteraction_MultiplyWithInteractionPotential(dOrb, pot, orb)
  hartree = Grid_InnerProduct(orb, dOrb)
  hartreeAnalytic = sqrt(2.0_R64) * PI**2.5_R64 * sigma**5

  call check(error, real(hartree, kind=R64), hartreeAnalytic, thr=1e-7_R64 * hartreeAnalytic)
  if (allocated(error)) error stop "T_SysInteraction_Prolate Hartree integral failure"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "all prolate interaction checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
