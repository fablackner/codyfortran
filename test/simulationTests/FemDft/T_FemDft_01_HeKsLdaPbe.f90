!> ctest-requires: libxc
!> He on grid.fem (one center, Galerkin point nucleus): Kohn-Sham LDA and PBE
!> ground states (groundSolver.scf.stdImpl with the Kohn-Sham mean field of
!> sysXc.fem, LOBPCG with the exact kinetic preconditioner) against the
!> converged radial-grid values (NIST LDA -2.834836; PBE -2.8929349 of
!> T_He3d_02 on a Ylm grid). The FEM nucleus is not variational and the
!> total energies carry its mesh error, which largely cancels in
!> E_PBE - E_LDA: that difference checks the FEM density gradients and the
!> GGA potential against the independent radial/angular gradients of
!> sysXc.ylm.
program T_FemDft_01_HeKsLdaPbe
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysXc
  use M_Orbs
  use M_Method_Mb_OrbBased_Tdks
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: E_LDA_REFERENCE = -2.834836_R64
  real(R64), parameter :: E_PBE_REFERENCE = -2.8929349_R64

  type(error_type), allocatable :: error
  real(R64) :: convThresh, epsThresh, energyLda, energyPbe, epsLda(1), epsPbe(1)
  integer(I32) :: maxScfIterations, nScfLda, nScfPbe, maxWarmLda, maxWarmPbe
  logical :: convergedLda, convergedPbe

  call Say_Hello
  call SetupFemScf("test/simulationTests/FemDft/T_FemDft_01_HeKsLdaPbe.json")
  convThresh = Json_Get("program.convThresh", 1e-11_R64)
  epsThresh = Json_Get("program.epsThresh", 1e-7_R64)
  maxScfIterations = Json_Get("program.maxScfIterations", 40)
  write (*, '(1X, A, I0)') "nodes ", Grid_nPoints

  energyLda = huge(1.0_R64)
  epsLda = huge(1.0_R64)
  call RunScf(nScfLda, maxWarmLda, convergedLda, energyLda, epsLda, "LDA", maxScfIterations, convThresh, epsThresh, &
              firstStepWarmQ=.false.)
  write (*, '(1X, A, 3F16.10)') "LDA: E_H, E_xc ", Method_Mb_OrbBased_Tdks_hartreeEnergy, Method_Mb_OrbBased_Tdks_xcEnergy

  ! PBE from the LDA orbitals and potential
  call jsonFile % remove("sysXc.fem.libxc.functionals")
  call jsonFile % add("sysXc.fem.libxc.functionals", ["gga_x_pbe", "gga_c_pbe"])
  call SysXc_Fabricate
  call SysXc_Setup
  energyPbe = energyLda
  epsPbe = epsLda
  call RunScf(nScfPbe, maxWarmPbe, convergedPbe, energyPbe, epsPbe, "PBE", maxScfIterations, convThresh, epsThresh, &
              firstStepWarmQ=.true.)

  print *
  print '(1X, A, F16.10, A, F14.8, A, ES10.2)', "He LDA (fem): E = ", energyLda, ", eps = ", epsLda(1), &
    ", E - reference ", energyLda - E_LDA_REFERENCE
  print '(1X, A, F16.10, A, F14.8, A, ES10.2)', "He PBE (fem): E = ", energyPbe, ", eps = ", epsPbe(1), &
    ", E - reference ", energyPbe - E_PBE_REFERENCE
  print '(1X, A, F14.10, A, F14.10, A, ES10.2)', "E_PBE - E_LDA: fem ", energyPbe - energyLda, ", radial ", &
    E_PBE_REFERENCE - E_LDA_REFERENCE, ", difference ", (energyPbe - energyLda) - (E_PBE_REFERENCE - E_LDA_REFERENCE)

  call check(error, convergedLda .and. convergedPbe)
  if (allocated(error)) error stop "T_FemDft_01_HeKsLdaPbe failure: a LOBPCG solve did not converge"
  call check(error, energyLda, E_LDA_REFERENCE, thr=1e-4_R64)
  if (allocated(error)) error stop "T_FemDft_01_HeKsLdaPbe failure: LDA energy"
  call check(error, energyPbe, E_PBE_REFERENCE, thr=1e-4_R64)
  if (allocated(error)) error stop "T_FemDft_01_HeKsLdaPbe failure: PBE energy"
  call check(error, energyPbe - energyLda, E_PBE_REFERENCE - E_LDA_REFERENCE, thr=5e-6_R64)
  if (allocated(error)) error stop "T_FemDft_01_HeKsLdaPbe failure: E_PBE - E_LDA differs from the radial grid"

  call Say_Goodbye

contains

  include "../Fem/FemScf.inc"

end program
