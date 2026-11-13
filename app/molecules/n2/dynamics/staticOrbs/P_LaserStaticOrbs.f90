!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> N2 in a laser pulse along the bond (sysGauge.prolate.length) in the fixed
!> basis of the MCSCF orbitals of ../../groundstate/mcscf: exact TDCI
!> (method tdci, full CI of the orbitals) or TD-2RDM (method
!> td2rdmStaticOrbs, reconstructed 3-RDM, optional purification). The JSON
!> file (first argument) selects the method; the programs share the
!> observables, so that the TD-2RDM runs compare directly with TDCI.
!>
!> Every program.outputStep the program writes to <program.outputName>.dat:
!>   t, E(t), dipole <z>, energy, natural occupations (spin up),
!>   min eigenvalue of D (up-down block), min eigenvalue of Q (up-down block)
!> The dipole is <z> = Tr(z rdm1) of the electrons (z from the bond center,
!> z_ij = <phi_i|z|phi_j> of the static orbitals), the energy the electronic
!> energy Tr(h1 rdm1) + Tr(h2 rdm2) with the field-free h1.
!>
!> program.kick (default 0): the initial state is replaced by
!> exp(-i kick Z) |Psi>, Z = sum z_ij a+_i a_j (TDCI: 12 Taylor terms;
!> TD-2RDM: its 2-RDM), a dipole kick for the linear response; it replaces
!> the pulse (no sysGauge block).
!>
!> TD-2RDM: with a purification block every time step ends with
!> Purification_PurifyTwoRdmState; with reconstruction.singlet.anchor the
!> MCSCF ground state is a stationary point of the TD-2RDM equations, and
!> anchor.nTransitions > 0 adds the spectral anchor (the transition RDMs of
!> the ground-state program with program.anchor).
!>
!> Inputs: P_LaserStaticOrbs.json (TDCI, Lanczos integrator: RK4 damps the CI
!> vector at |E| = 133 Ha), td2rdmNy.json (NY with ccQ), td2rdmNyAnchor.json
!> (anchored), td2rdmNySpectral.json (spectral anchor, D/Q purification).
!> Results in app/README.md.
program P_LaserStaticOrbs
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_RdmDiagonalize
  use M_Utils_RdmHole
  use M_Utils_RdmObservables
  use M_Utils_LapackLib
  use M_Grid
  use M_Grid_Prolate
  use M_SysKinetic
  use M_SysPotential
  use M_SysGauge
  use M_SysGauge_Length
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_TwoRdm
  use M_TwoRdmInit
  use M_Reconstruction
  use M_Purification
  use M_Method
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator

  implicit none

  real(R64)    :: time, endTime, outputStep, timeStep, kick
  integer(I32) :: nTimeSteps, innerStep, io
  logical      :: td2rdmQ, purifyQ, gaugeQ
  character(len=:), allocatable :: outputName
  complex(R64), allocatable :: zMat(:, :), h1Free(:, :), h2Free(:, :, :, :)
  real(R64), allocatable :: z(:)
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile()

  td2rdmQ = Json_GetExistence("method.mb.orbBased.td2rdmStaticOrbs")
  purifyQ = Json_GetExistence("purification")
  gaugeQ = Json_GetExistence("sysGauge")

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  if (gaugeQ) call SysGauge_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  if (td2rdmQ) then
    call TwoRdm_Fabricate
    call TwoRdmInit_Fabricate
    call Reconstruction_Fabricate
    if (purifyQ) call Purification_Fabricate
  else
    call ConfigList_Fabricate
    call Coeffs_Fabricate
    call CoeffsInit_Fabricate
  end if
  IntegratorListInput(1) % TimeDerivative => Method_TimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Say_Fabricate("program")
  endTime = Json_Get("program.endTime", 100.0_R64)
  outputStep = Json_Get("program.outputStep", 0.2_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 20)
  outputName = Json_Get("program.outputName", "observables")
  kick = Json_Get("program.kick", 0.0_R64)
  if (abs(kick) > 0.0_R64 .and. gaugeQ) error stop "program.kick replaces the pulse: remove the sysGauge block"

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  if (gaugeQ) call SysGauge_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  if (td2rdmQ) then
    call TwoRdm_Setup
    call TwoRdmInit_Setup
  else
    call ConfigList_Setup
    call Coeffs_Setup
    call CoeffsInit_Setup
  end if
  call Method_Setup
  if (td2rdmQ) then
    call Reconstruction_Setup
    if (purifyQ) call Purification_Setup
  end if
  call IntegratorList_Setup
  call Propagator_Setup

  ! field-free h1 and h2, dipole matrix of the static orbitals
  call Method_Mb_OrbBased_FillH0(h1Free, Orbs_orbs, 0.0_R64)
  call Method_Mb_OrbBased_FillH2(h2Free, Orbs_orbs, 0.0_R64)
  call Grid_Prolate_FillNodeCoordinate(z, [0.0_R64, 0.0_R64, 1.0_R64])
  call Method_Mb_OrbBased_FillOneBodyMatrix(zMat, Orbs_orbs, spread(z, 2, Orbs_nOrbsInState) * Orbs_orbs)

  if (abs(kick) > 0.0_R64) then
    if (td2rdmQ) then
      call KickTwoRdm(kick)
    else
      call KickCoeffs(kick)
    end if
  end if

  !==========================================
  call Say_Section("propagation")
  !==========================================

  open (newunit=io, file=outputName//".dat", status="replace")

  time = 0.0_R64
  timeStep = outputStep / nTimeSteps
  call WriteObservables(time)

  do while (time < endTime - 0.5_R64 * timeStep)
    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
      if (td2rdmQ .and. purifyQ) call Purification_PurifyTwoRdmState(TwoRdm_twoRdmState, h2Free)
    end do
    call WriteObservables(time)
  end do

  close (io)

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> rdm2 of exp(-i kick Z)|Psi>: the 2-RDM is ket-first,
  !> rdm2(i,j,k,l) = <a+_k a+_l a_j a_i>, and exp(i kick Z) a_i exp(-i kick Z)
  !> = sum_p M(i,p) a_p with M = exp(-i kick z) (z the spin-orbital dipole
  !> matrix), so rdm2(i,j,k,l) <- sum M(i,p) M(j,q) rdm2(p,q,r,s)
  !> conjg(M(k,r) M(l,s)).
  subroutine KickTwoRdm(kick)
    use M_Utils_RdmBasisChange

    real(R64), intent(in) :: kick

    complex(R64), allocatable :: rdm2(:, :, :, :), m(:, :), evecs(:, :)
    real(R64), allocatable :: evals(:)
    integer(I32) :: k

    call LapackLib_DiagonalizeGeneric(evals, evecs, 0.5_R64 * (zMat + conjg(transpose(zMat))), .true.)
    allocate (m, mold=evecs)
    do k = 1, size(evals)
      m(:, k) = evecs(:, k) * exp(cmplx(0.0_R64, -kick * evals(k), R64))
    end do
    m = matmul(m, conjg(transpose(evecs)))

    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, TwoRdm_twoRdmState)
    call RdmBasisChange_Transform(rdm2, 4, m)
    call TwoRdm_FillTwoRdmStateFromRdm2(TwoRdm_twoRdmState, rdm2)

  end subroutine

  !> Coeffs_coeffs <- exp(-i kick Z) Coeffs_coeffs by 12 Taylor terms
  !> (kick |Z| << 1).
  subroutine KickCoeffs(kick)
    real(R64), intent(in) :: kick

    complex(R64), allocatable :: term(:), zTerm(:)
    integer(I32) :: k

    term = Coeffs_coeffs
    allocate (zTerm, mold=term)
    do k = 1, 12
      zTerm = 0.0_R64
      call Coeffs_ApplyH1FillRdm1(term, dCoeffs_=zTerm, h1_=zMat)
      term = cmplx(0.0_R64, -kick / k, R64) * zTerm
      Coeffs_coeffs = Coeffs_coeffs + term
    end do

  end subroutine

  subroutine WriteObservables(t)
    real(R64), intent(in) :: t

    complex(R64), allocatable :: rdm1(:, :), rdm2(:, :, :, :), natOrbs(:, :)
    real(R64), allocatable :: natOcc(:)
    real(R64) :: energy, dipole, minD, minQ
    integer(I32) :: n

    call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
    call Method_Mb_OrbBased_FillRdm2(rdm2, Method_state)

    energy = RdmObservables_Energy(rdm1, h1Free, rdm2, h2Free)
    dipole = real(sum(zMat * transpose(rdm1)), R64)

    n = Orbs_nOrbsInState
    call RdmDiagonalize_Rdm1(natOcc, natOrbs, rdm1(1:n, 1:n))
    call RdmHole_MinEigenvaluesDQab(minD, minQ, rdm2(1:n, n + 1:, 1:n, n + 1:))

    write (io, '(F10.4, ES16.8, F20.12, F24.14, *(ES16.8))') t, SysGauge_Length_FieldAmplitude(t), &
      dipole, energy, natOcc, minD, minQ
    flush (io)
    write (*, '(A,F9.3,A,F14.8,A,F20.12,A,2ES11.3)') " t =", t, "  dipole =", dipole, "  energy =", energy, &
      "  min eig D, Q =", minD, minQ

  end subroutine

end program
