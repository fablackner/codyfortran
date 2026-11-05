!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Restricted Hartree–Fock ground state of a molecule on the finite-element
!> grid (grid.fem). Run commands in app/README.md. The JSON file is the first
!> command-line argument; every file is written to the working directory.
!>
!> Wiring (as test/simulationTests/Fem/T_Fem_03 and T_Fem_04): method tdhx
!> with restricted orbitals (orbs.restrictedQ, nOrbs = nBodies = number of
!> doubly occupied orbitals per spin), groundSolver.scf.stdImpl (Hartree
!> potential of the total density, adaptively compressed exchange with aceQ,
!> mixing of the Hartree potential) and the eigensolver diagonalizerList.lobpcg
!> in the FEM metric (metricWeights = Grid_Fem_weights), preconditioned by the
!> exact FDM solve (T + max(−θ, σ_min))⁻¹ (GroundSolver_PreconditionAction).
!> The start orbitals of orbsInit (orbsInit.fem.sto, or orbsInit.load to
!> continue an earlier run on the same mesh) are orthonormalized first.
!>
!> SCF loop: one GroundSolver_Approach per step (Fock operator of the current
!> orbitals, its lowest eigenvectors, mixing) until |ΔE_el| <
!> program.convThresh and, if program.epsThresh > 0, max|Δε| <
!> program.epsThresh (orbital energies converge only like √ΔE), within
!> program.maxIter steps. One line per step: E_el, ΔE, max|Δε|, LOBPCG
!> iterations and Fock products, wall time of the step.
!>
!> Results: electronic energy (Method_GetEnergy), nuclear repulsion
!> (SysPotential_Fem_Coulomb_nuclearRepulsion: point charges
!> sysPotential.fem.coulomb at the grid.fem centers), total energy, occupied
!> orbital energies, electron number Tr(ρS) and the dipole moment
!> μ = Σ_I Z_I R_I − ⟨r⟩_el (nuclear + electronic part; independent of the
!> origin for a neutral molecule).
!>
!> Files:
!>   orb01_01.in ...  the orbitals (Orbs_SaveOrbs): initial state of
!>                    P_KickTdhfFem (orbsInit.load) on the IDENTICAL mesh
!>   femMesh.dat      nLocals and element bounds of the mesh; P_KickTdhfFem
!>                    stops if its own mesh differs
!> If the SCF does not converge within program.maxIter, the program says so,
!> still writes both files (continue with orbsInit.load) and stops with an
!> error code.
program P_GroundScfFem
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_PrinterObservableFem
  use M_Grid
  use M_Grid_Fem
  use M_SysKinetic
  use M_SysPotential
  use M_SysPotential_Fem_Coulomb, only: SysPotential_Fem_Coulomb_nuclearRepulsion, SysPotential_Fem_Coulomb_charges
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb
  use M_Method_Mb_OrbBased
  use M_DiagonalizerList
  use M_DiagonalizerList_Lobpcg, only: T_DiagonalizerList_E_Lobpcg
  use M_Mixing
  use M_GroundSolver
  use M_GroundSolver_Scf

  implicit none

  real(R64), parameter :: HARTREE_EV = 27.211386245988_R64
  real(R64), parameter :: AU_DEBYE = 2.541746473_R64

  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)
  complex(R64), allocatable :: rdm1(:, :)
  real(R64), allocatable :: eps(:), epsOld(:)
  real(R64) :: convThresh, epsThresh, energy, energyOld, secondsStep, secondsScf
  real(R64) :: nElectrons, electronPosition(3), nuclearDipole(3), dipole(3)
  integer(I64) :: count0, countStep
  integer(I32) :: maxIter, nOcc, iScf, nScf, nIter, nMatVec, nMatVecTotal, i
  logical :: convergedQ
  character(len=10) :: dEnergyText, dEpsText

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile()

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  call Mixing_Fabricate
  call GroundSolver_Fabricate
  DiagonalizerListInput(1) % ApplyMatOnVec => GroundSolver_Scf_HartreeFockAction
  DiagonalizerListInput(1) % ApplyPrecOnVec => GroundSolver_PreconditionAction
  DiagonalizerListInput(1) % dim = Grid_nPoints
  call DiagonalizerList_Fabricate(DiagonalizerListInput)

  call Say_Fabricate("program")
  convThresh = Json_Get("program.convThresh", 1e-8_R64)
  epsThresh = Json_Get("program.epsThresh", 0.0_R64)
  maxIter = Json_Get("program.maxIter", 50)
  if (.not. (convThresh > 0.0_R64)) error stop "P_GroundScfFem: program.convThresh must be positive"
  if (maxIter < 2) error stop "P_GroundScfFem: program.maxIter must be at least 2 (|dE| needs two SCF steps)"

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call GroundSolver_Setup
  call DiagonalizerList_Setup

  ! LOBPCG works in the FEM metric (diagonal mass, BLAS-3 Gram matrices)
  DiagonalizerList(1) % e % metricWeights = Grid_Fem_weights

  ! OrbsInit only normalizes: orthonormalize the start orbitals
  call Orbs_Orthonormalize(Orbs_orbs)

  ! Doubly occupied spatial orbitals (restricted closed shell)
  nOcc = GroundSolver_NumberOfSpinUpOrbs()
  if (DiagonalizerList(1) % e % nEvals < nOcc) &
    error stop "P_GroundScfFem: the diagonalizer must provide at least one eigenvalue per occupied orbital (nEvals)"
  allocate (eps(nOcc), epsOld(nOcc))

  !==========================================
  call Say_Section("self-consistent field")
  !==========================================

  call PrintSystem

  call system_clock(count0)
  energy = 0.0_R64
  eps(:) = 0.0_R64
  nMatVecTotal = 0
  convergedQ = .false.
  nScf = 0

  do iScf = 1, maxIter
    nScf = iScf
    energyOld = energy
    epsOld(:) = eps(:)

    call system_clock(countStep)
    call GroundSolver_Approach(Method_state, time=0.0_R64)
    energy = Method_GetEnergy(0.0_R64)
    eps(:) = DiagonalizerList(1) % e % evals(1:nOcc)
    secondsStep = SecondsSince(countStep)

    nIter = -1
    nMatVec = -1
    select type (e => DiagonalizerList(1) % e)
    type is (T_DiagonalizerList_E_Lobpcg)
      nIter = e % nIter
      nMatVec = e % nMatVec
      nMatVecTotal = nMatVecTotal + nMatVec
    end select

    dEnergyText = "         -"
    dEpsText = "         -"
    if (iScf > 1) then
      write (dEnergyText, '(ES10.2)') energy - energyOld
      write (dEpsText, '(ES10.2)') maxval(abs(eps - epsOld))
    end if
    write (*, '(1X, A, I3, A, F20.12, 4A, A, I4, A, I5, A, L1, A, F9.2, A)') "scf ", iScf, ": E_el ", energy, &
      "  dE ", dEnergyText, "  max|deps| ", dEpsText, "  lobpcg iterations ", nIter, "  Fock products ", nMatVec, &
      "  converged ", DiagonalizerList(1) % e % convergedQ, "  time ", secondsStep, " s"
    flush (6)

    if (iScf > 1) then
      if (abs(energy - energyOld) < convThresh .and. &
          (epsThresh <= 0.0_R64 .or. maxval(abs(eps - epsOld)) < epsThresh)) then
        convergedQ = .true.
        exit
      end if
    end if
  end do

  secondsScf = SecondsSince(count0)
  write (*, *)
  write (*, '(1X, A, I0, A, I0, A, F10.2, A, F8.2, A)') "SCF: ", nScf, " steps, ", nMatVecTotal, &
    " Fock matrix-vector products, ", secondsScf, " s (", secondsScf / nScf, " s per step)"
  if (.not. convergedQ) then
    write (*, '(1X, A, I0, A, ES10.2, A, ES10.2, A)') "WARNING: SCF NOT CONVERGED within program.maxIter = ", maxIter, &
      " steps (last |dE| = ", abs(energy - energyOld), ", program.convThresh = ", convThresh, ")"
  end if

  !==========================================
  call Say_Section("printout results")
  !==========================================

  call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
  call PrinterObservableFem_DumpNorm(rdm1, Orbs_orbs, "", .false., 0.0_R64, norm_=nElectrons)
  call PrinterObservableFem_DumpDipole(rdm1, Orbs_orbs, "", .false., 0.0_R64, dipole_=electronPosition)
  nuclearDipole(:) = matmul(Grid_Fem_centers, SysPotential_Fem_Coulomb_charges)
  dipole(:) = nuclearDipole(:) - electronPosition(:)

  write (*, *)
  write (*, '(1X, A, F22.12, A)') "electronic energy   ", energy, " Ha"
  write (*, '(1X, A, F22.12, A)') "nuclear repulsion   ", SysPotential_Fem_Coulomb_nuclearRepulsion, " Ha"
  write (*, '(1X, A, F22.12, A)') "total energy        ", energy + SysPotential_Fem_Coulomb_nuclearRepulsion, " Ha"
  write (*, *)
  write (*, '(1X, A)') "occupied orbital energies (Ha, eV):"
  do i = 1, nOcc
    write (*, '(1X, A, I3, F20.10, F16.6)') "  ", i, eps(i), eps(i) * HARTREE_EV
  end do
  write (*, *)
  write (*, '(1X, A, F18.12, A, I0, A)') "electrons Tr(rho S) ", nElectrons, " (", Method_Mb_nBodiesSum, " bodies)"
  write (*, *)
  write (*, '(1X, A)') "dipole moment mu = sum_I Z_I R_I - <r>_el (a.u.):         x               y               z"
  write (*, '(1X, A, 3F16.10)') "  nuclear     sum_I Z_I R_I               ", nuclearDipole
  write (*, '(1X, A, 3F16.10)') "  electronic  -<r>_el                     ", -electronPosition
  write (*, '(1X, A, 3F16.10)') "  total       mu                          ", dipole
  write (*, '(1X, A, F14.10, A, F14.10, A)') "  |mu| = ", norm2(dipole), " a.u. = ", norm2(dipole) * AU_DEBYE, " D"
  if (abs(sum(SysPotential_Fem_Coulomb_charges) - nElectrons) > 1e-6_R64) &
    write (*, '(1X, A)') "  (charged system: the dipole moment depends on the origin)"

  !==========================================
  call Say_Section("save data")
  !==========================================

  call Orbs_SaveOrbs(Orbs_orbs)
  call WriteFemMeshFile("femMesh.dat")
  write (*, '(1X, A, I0, A)') "saved ", Orbs_nOrbsInState, " orbitals (orbXX_YY.in) and femMesh.dat"

  if (.not. convergedQ) &
    error stop "P_GroundScfFem: SCF not converged within program.maxIter (orbitals of the last step saved; continue with orbsInit.load)"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Mesh size, nuclei and electrons of the run.
  subroutine PrintSystem
    integer(I32) :: iCenter, iAxis

    write (*, '(1X, A, I0, A, 3(I0, 1X), A, I0, A, 3(I0, 1X), A)') "mesh: nLocals ", Grid_Fem_nLocals, &
      ", elements ", (Grid_Fem_axes(iAxis) % nElements, iAxis=1, 3), ", interior nodes N = ", Grid_nPoints, &
      " = ", Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, "(x y z)"
    write (*, '(1X, A)') "nuclei (charge, position):"
    do iCenter = 1, Grid_Fem_nCenters
      write (*, '(1X, A, F8.3, 3F16.10)') "  ", SysPotential_Fem_Coulomb_charges(iCenter), Grid_Fem_centers(:, iCenter)
    end do
    write (*, '(1X, A, I0, A, I0, A)') "electrons: ", Method_Mb_nBodiesSum, " in ", nOcc, &
      " doubly occupied orbitals"
    write (*, *)

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Writes nLocals and the element bounds of the three axes (full precision),
  !> the fingerprint of the mesh the saved orbitals live on.
  subroutine WriteFemMeshFile(filename)
    character(len=*), intent(in) :: filename

    integer(I32) :: io, iAxis

    open (newunit=io, file=filename, status="replace", action="write")
    write (io, '(A)') "# grid.fem mesh of the orbitals orbXX_YY.in (P_GroundScfFem); checked by P_KickTdhfFem"
    write (io, '(A)') "# nLocals, then per axis x, y, z: nElements and the element bounds"
    write (io, '(I0)') Grid_Fem_nLocals
    do iAxis = 1, 3
      write (io, '(I0)') Grid_Fem_axes(iAxis) % nElements
      write (io, '(*(ES25.17))') Grid_Fem_axes(iAxis) % bounds
    end do
    close (io)

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Wall-clock seconds since count (system_clock counts).
  function SecondsSince(count) result(res)
    integer(I64), intent(in) :: count
    real(R64) :: res

    integer(I64) :: countNow, rate

    call system_clock(countNow, rate)
    res = real(countNow - count, kind=R64) / real(rate, kind=R64)

  end function

end program
