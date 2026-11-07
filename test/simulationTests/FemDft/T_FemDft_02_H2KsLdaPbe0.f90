!> ctest-requires: libxc
!> H2 (R = 1.4) on the grid.fem mesh of T_Fem_03: Kohn-Sham ground states
!> with LDA (Slater + VWN5) and with the global hybrid PBE0 (a = 1/4 exact
!> exchange of SysInteraction, applied by ACE and by the exact exchange: both
!> must give the same energy and orbital energy); writes the LDA orbital
!> orb01_01.in for the TDDFT test T_FemDft_03. DIIS on the orbitals
!> (mixing.diis without start gate, startThreshold 0: the complete iterate also for the hybrid,
!> whose exact exchange would lag a mixed potential), reset between the
!> three SCFs (Mixing_Reset: each has its own fixed-point map).
program T_FemDft_02_H2KsLdaPbe0
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysXc
  use M_SysPotential_Fem_Coulomb, only: SysPotential_Fem_Coulomb_nuclearRepulsion
  use M_Orbs
  use M_Mixing
  use M_GroundSolver
  use M_Method_Mb_OrbBased_Tdks
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  real(R64) :: convThresh, epsThresh, energyLda, energyAce, energyExact, epsLda(1), epsAce(1), epsExact(1)
  integer(I32) :: maxScfIterations, nScf, maxWarm
  logical :: convergedQ, allConvergedQ, foundQ

  call Say_Hello
  call SetupFemScf("test/simulationTests/FemDft/T_FemDft_02_H2KsLdaPbe0.json")
  convThresh = Json_Get("program.convThresh", 1e-12_R64)
  epsThresh = Json_Get("program.epsThresh", 1e-8_R64)
  maxScfIterations = Json_Get("program.maxScfIterations", 60)
  write (*, '(1X, A, I0)') "nodes ", Grid_nPoints

  ! LDA
  energyLda = huge(1.0_R64)
  epsLda = huge(1.0_R64)
  call RunScf(nScf, maxWarm, convergedQ, energyLda, epsLda, "LDA", maxScfIterations, convThresh, epsThresh, &
              firstStepWarmQ=.false.)
  allConvergedQ = convergedQ
  call Orbs_SaveOrbs(Orbs_orbs)

  ! PBE0 with ACE
  call jsonFile % remove("sysXc.fem.libxc.functionals")
  call jsonFile % add("sysXc.fem.libxc.functionals", ["hyb_gga_xc_pbeh"])
  call SysXc_Fabricate
  call SysXc_Setup
  call Mixing_Reset
  energyAce = energyLda
  epsAce = epsLda
  call RunScf(nScf, maxWarm, convergedQ, energyAce, epsAce, "PBE0, ACE", maxScfIterations, convThresh, epsThresh, &
              firstStepWarmQ=.true.)
  allConvergedQ = allConvergedQ .and. convergedQ
  write (*, '(1X, A, 4F16.10)') "PBE0: E_H, E_xc, E_x, a ", Method_Mb_OrbBased_Tdks_hartreeEnergy, &
    Method_Mb_OrbBased_Tdks_xcEnergy, Method_Mb_OrbBased_Tdks_exactExchangeEnergy, SysXc_exactExchangeFraction

  ! PBE0 with the exact exchange
  call jsonFile % update("groundSolver.scf.stdImpl.aceQ", .false., foundQ)
  if (.not. foundQ) error stop "T_FemDft_02_H2KsLdaPbe0: groundSolver.scf.stdImpl.aceQ not found"
  call GroundSolver_Fabricate
  call Mixing_Reset
  energyExact = energyAce
  epsExact = epsAce
  call RunScf(nScf, maxWarm, convergedQ, energyExact, epsExact, "PBE0, exact exchange", maxScfIterations, &
              convThresh, epsThresh, firstStepWarmQ=.true.)

  print *
  print '(1X, A, F18.12, A, F18.12, A, F14.10)', "H2 LDA:  E_el ", energyLda, ", E_tot ", &
    energyLda + SysPotential_Fem_Coulomb_nuclearRepulsion, ", eps ", epsLda(1)
  print '(1X, A, F18.12, A, F18.12, A, F14.10)', "H2 PBE0: E_el ", energyAce, ", E_tot ", &
    energyAce + SysPotential_Fem_Coulomb_nuclearRepulsion, ", eps ", epsAce(1)
  print '(1X, A, 2ES12.3)', "PBE0 exact exchange - ACE: dE, deps ", energyExact - energyAce, epsExact(1) - epsAce(1)

  call check(error, allConvergedQ)
  if (allocated(error)) error stop "T_FemDft_02_H2KsLdaPbe0 failure: a LOBPCG solve with ACE did not converge"
  call check(error, energyLda, -1.851804662_R64, thr=1e-8_R64)
  if (allocated(error)) error stop "T_FemDft_02_H2KsLdaPbe0 failure: LDA energy (regression)"
  call check(error, energyExact, energyAce, thr=1e-9_R64)
  if (allocated(error)) error stop "T_FemDft_02_H2KsLdaPbe0 failure: PBE0 energy, exact exchange vs ACE"
  call check(error, epsExact(1), epsAce(1), thr=1e-7_R64)
  if (allocated(error)) error stop "T_FemDft_02_H2KsLdaPbe0 failure: PBE0 orbital energy, exact exchange vs ACE"

  call Say_Goodbye

contains

  include "../Fem/FemScf.inc"

end program
