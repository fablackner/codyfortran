!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Exact diagonalization of the TDCI Hamiltonian of a Hubbard lattice, shared
! by the LAPACK and ARPACK tests of FermiHubbard (coeffs.generic and
! coeffs.hubbard) and BoseHubbard; included after `contains`.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Loads jsonFileName, computes the lowest eigenpair of
!> Method_ApplyHamiltonian with the first diagonalizer of the input and stops
!> unless the energy of the eigenvector (Method_GetEnergy) agrees with the
!> eigenvalue and with referenceEnergy to 1e-7.
subroutine CheckExactDiagGroundState(jsonFileName, referenceEnergy)
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Printer_Spectrum
  use M_Grid
  use M_SysKinetic
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_DiagonalizerList
  use testdrive, only: check, error_type

  character(len=*), intent(in) :: jsonFileName
  real(R64), intent(in) :: referenceEnergy

  real(R64) :: energy
  type(error_type), allocatable :: error
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)
  logical :: configListQ

  !==========================================
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)
  ! coeffs.hubbard has no configuration list
  configListQ = Json_GetExistence("configList")

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  if (configListQ) call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  DiagonalizerListInput(1) % ApplyMatOnVec => Method_ApplyHamiltonian
  DiagonalizerListInput(1) % dim = Coeffs_nCoeffs
  call DiagonalizerList_Fabricate(DiagonalizerListInput)

  call Grid_Setup
  call SysKinetic_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  if (configListQ) call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call DiagonalizerList_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  call DiagonalizerList(1) % e % Diagonalize(0.0_R64, .true.)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  if (DiagonalizerList(1) % e % nFound .eq. 0) error stop "CheckExactDiagGroundState: no eigenvalue found"

  Method_state(:) = DiagonalizerList(1) % e % evecs(:, 1)
  energy = Method_GetEnergy(0.0_R64)

  write (*, *)
  write (*, *) "groundstate energy: ", energy, DiagonalizerList(1) % e % evals(1)
  write (*, *)

  call Printer_Spectrum_DumpSpectrum("eigenspectrum.dat", .true.)

  call check(error, energy, DiagonalizerList(1) % e % evals(1), thr=1e-7_R64)
  if (allocated(error)) error stop "CheckExactDiagGroundState: energy of the eigenvector differs from the eigenvalue"

  call check(error, energy, referenceEnergy, thr=1e-7_R64)
  if (allocated(error)) error stop "CheckExactDiagGroundState: ground-state energy off the reference"

end subroutine
