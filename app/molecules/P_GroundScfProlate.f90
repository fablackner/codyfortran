!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Restricted Hartree–Fock ground state of a homonuclear diatomic molecule on
!> the prolate spheroidal grid (grid.prolate, nuclei at the foci z = ±R/2).
!> Run commands in app/README.md. The JSON file is the first command-line
!> argument; every file is written to the working directory.
!>
!> Wiring: method tdhx with restricted orbitals (orbs.restrictedQ, nOrbs =
!> nBodies = number of doubly occupied orbitals per spin),
!> groundSolver.scf.stdImpl (Fock operator of the current orbitals, DIIS on
!> the orbitals: mixTarget "orbitals", mixing.diis with startThreshold 0) and
!> one diagonalizer for the Fock operator on the whole grid (all azimuthal
!> channels m = 0, ±1, ... up to grid.prolate.mmax; σ and π orbitals come out
!> of the same eigenproblem). A lobpcg diagonalizer works in the prolate
!> metric (metricWeights = Grid_weights), preconditioned by the exact
!> (T + max(−θ, σ_min))⁻¹ of the prolate Laplacian
!> (GroundSolver_PreconditionAction); use it with aceQ for heavier nuclei
!> (graded xi elements), ARPACK for H2. The start
!> orbitals of orbsInit (orbsInit.prolate.lcao) are orthonormalized first.
!>
!> SCF loop: one GroundSolver_Approach per step until |ΔE_el| <
!> program.convThresh and, if program.epsThresh > 0, max|Δε| <
!> program.epsThresh, within program.maxIter steps.
!>
!> Results: electronic energy, nuclear repulsion Z²/R, total energy, and per
!> occupied orbital its energy and symmetry label (|m| of the dominant
!> azimuthal channel: σ, π, δ; inversion parity g/u from
!> (ξ, η, φ) → (ξ, −η, φ + π)).
!>
!> Files: orb01_01.in ... (Orbs_SaveOrbs). If the SCF does not converge within
!> program.maxIter, the program says so, still writes the orbitals (continue
!> with orbsInit.load) and stops with an error code.
program P_GroundScfProlate
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Prolate
  use M_SysKinetic
  use M_SysPotential
  use M_SysPotential_Prolate_Coulomb, only: SysPotential_Prolate_Coulomb_charge
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb
  use M_DiagonalizerList
  use M_DiagonalizerList_Lobpcg, only: T_DiagonalizerList_E_Lobpcg
  use M_Mixing
  use M_GroundSolver
  use M_GroundSolver_Scf

  implicit none

  real(R64), parameter :: HARTREE_EV = 27.211386245988_R64

  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)
  real(R64), allocatable :: eps(:), epsOld(:)
  real(R64) :: convThresh, epsThresh, energy, energyOld, nuclearRepulsion, secondsStep, secondsScf
  integer(I64) :: count0, countStep
  integer(I32) :: maxIter, nOcc, iScf, nScf, nMatVec, nMatVecTotal, i
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
  convThresh = Json_Get("program.convThresh", 1e-10_R64)
  epsThresh = Json_Get("program.epsThresh", 0.0_R64)
  maxIter = Json_Get("program.maxIter", 50)
  if (.not. (convThresh > 0.0_R64)) error stop "P_GroundScfProlate: program.convThresh must be positive"
  if (maxIter < 2) error stop "P_GroundScfProlate: program.maxIter must be at least 2 (|dE| needs two SCF steps)"

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

  ! LOBPCG works in the prolate metric (quadrature x Jacobian); LAPACK and
  ! ARPACK ignore it
  DiagonalizerList(1) % e % metricWeights = Grid_weights

  ! OrbsInit only normalizes: orthonormalize the start orbitals
  call Orbs_Orthonormalize(Orbs_orbs)

  nOcc = GroundSolver_NumberOfSpinUpOrbs()
  if (DiagonalizerList(1) % e % nEvals < nOcc) &
    error stop "P_GroundScfProlate: the diagonalizer must provide at least one eigenvalue per occupied orbital (nEvals)"
  allocate (eps(nOcc), epsOld(nOcc))

  nuclearRepulsion = SysPotential_Prolate_Coulomb_charge**2 / Grid_Prolate_R

  !==========================================
  call Say_Section("self-consistent field")
  !==========================================

  write (*, '(1X, A, F10.6, A, F8.3, A, I0, A, I0, A, I0, A, I0, A)') "R = ", Grid_Prolate_R, ", Z = ", &
    SysPotential_Prolate_Coulomb_charge, ", grid points ", Grid_nPoints, " = ", Grid_Prolate_nXi, " (xi) x ", &
    Grid_Prolate_nEta, " (eta) x ", 2 * Grid_Prolate_mmax + 1, " (m)"
  write (*, '(1X, A, I0, A, I0, A)') "electrons: ", Method_Mb_nBodiesSum, " in ", nOcc, " doubly occupied orbitals"
  write (*, *)

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

    nMatVec = -1
    select type (e => DiagonalizerList(1) % e)
    type is (T_DiagonalizerList_E_Lobpcg)
      nMatVec = e % nMatVec
      nMatVecTotal = nMatVecTotal + nMatVec
    end select

    dEnergyText = "         -"
    dEpsText = "         -"
    if (iScf > 1) then
      write (dEnergyText, '(ES10.2)') energy - energyOld
      write (dEpsText, '(ES10.2)') maxval(abs(eps - epsOld))
    end if
    write (*, '(1X, A, I3, A, F20.12, 4A, A, I5, A, L1, A, F8.3, A)') "scf ", iScf, ": E_el ", energy, &
      "  dE ", dEnergyText, "  max|deps| ", dEpsText, "  Fock products ", nMatVec, &
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
  write (*, '(1X, A, I0, A, F10.3, A)') "SCF: ", nScf, " steps, ", secondsScf, " s"
  if (nMatVecTotal > 0) write (*, '(1X, A, I0)') "Fock matrix-vector products (lobpcg): ", nMatVecTotal
  if (.not. convergedQ) then
    write (*, '(1X, A, I0, A, ES10.2, A, ES10.2, A)') "WARNING: SCF NOT CONVERGED within program.maxIter = ", maxIter, &
      " steps (last |dE| = ", abs(energy - energyOld), ", program.convThresh = ", convThresh, ")"
  end if

  !==========================================
  call Say_Section("printout results")
  !==========================================

  write (*, *)
  write (*, '(1X, A, F22.12, A)') "electronic energy   ", energy, " Ha"
  write (*, '(1X, A, F22.12, A)') "nuclear repulsion   ", nuclearRepulsion, " Ha"
  write (*, '(1X, A, F22.12, A)') "total energy        ", energy + nuclearRepulsion, " Ha"
  write (*, *)
  write (*, '(1X, A)') "occupied orbital energies (Ha, eV), symmetry of the Fock eigenvectors:"
  do i = 1, nOcc
    write (*, '(1X, A, I3, F20.10, F16.6, 3X, A)') "  ", i, eps(i), eps(i) * HARTREE_EV, &
      SymmetryLabel(DiagonalizerList(1) % e % evecs(:, i))
  end do

  !==========================================
  call Say_Section("save data")
  !==========================================

  call Orbs_SaveOrbs(Orbs_orbs)
  write (*, '(1X, A, I0, A)') "saved ", Orbs_nOrbsInState, " orbitals (orbXX_YY.in)"

  if (.not. convergedQ) &
    error stop "P_GroundScfProlate: SCF not converged within program.maxIter (orbitals of the last step saved)"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Symmetry label of an orbital: |m| of the azimuthal channel with the
  !> largest weight (σ, π, δ, ...) and the inversion parity g/u. Inversion
  !> maps (ξ, η, φ) to (ξ, −η, φ + π): channel m picks up (−1)^m, and the
  !> Gauss–Legendre η nodes are symmetric (index j ↔ nEta + 1 − j).
  function SymmetryLabel(orb) result(label)
    complex(R64), intent(in) :: orb(:)
    character(len=8) :: label

    character(len=5), parameter :: NAMES(0:3) = ["sigma", "pi   ", "delta", "phi  "]
    real(R64) :: weightM, weightBest
    integer(I32) :: mChannel, mBest, mAbs, iXi, jEta, base, iPoint, iMirror
    complex(R64) :: overlap

    weightBest = -1.0_R64
    mBest = 0
    do mChannel = 0, 2 * Grid_Prolate_mmax
      base = mChannel * Grid_Prolate_nSpatial
      weightM = sum(Grid_weights(base + 1:base + Grid_Prolate_nSpatial) * &
                    abs(orb(base + 1:base + Grid_Prolate_nSpatial))**2)
      if (weightM > weightBest) then
        weightBest = weightM
        mBest = mChannel
      end if
    end do

    ! channel order m = 0, -1, +1, -2, +2, ...
    mAbs = (mBest + 1) / 2

    base = mBest * Grid_Prolate_nSpatial
    overlap = (0.0_R64, 0.0_R64)
    do jEta = 1, Grid_Prolate_nEta
      do iXi = 1, Grid_Prolate_nXi
        iPoint = base + (jEta - 1) * Grid_Prolate_nXi + iXi
        iMirror = base + (Grid_Prolate_nEta - jEta) * Grid_Prolate_nXi + iXi
        overlap = overlap + Grid_weights(iPoint) * conjg(orb(iPoint)) * orb(iMirror)
      end do
    end do

    label = "m"
    if (mAbs <= 3) label = NAMES(mAbs)
    if (real(overlap, kind=R64) * (-1.0_R64)**mAbs >= 0.0_R64) then
      label = trim(label)//"_g"
    else
      label = trim(label)//"_u"
    end if

  end function

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
