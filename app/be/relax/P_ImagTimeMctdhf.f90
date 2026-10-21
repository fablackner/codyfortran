!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> MCTDHF ground state of an atom (Ylm grid; Be with m = 0 orbitals as in
!> Lackner et al., PRA 95, 033414 (2017), or Ne, app/ne/hhg) by imaginary-time
!> propagation.
!>
!> Imaginary-time Strang splitting exp(-h/2) exp(-W) exp(-h/2):
!> - one-body part: every orbital is propagated with exp(-h dtau/2) at fixed
!>   CI coefficients, which is the exact one-body evolution of the many-body
!>   state. Without method.mb.orbBased.exactFieldFreeQ the Short Iterative
!>   Lanczos does this in imaginary time; with it, exp(-h0 dtau/4) of the
!>   field-free h0 is applied exactly on both sides of the SIL step, which
!>   then only sees the vanishing gauge term (ImagTimeLin);
!> - interaction part: classical RK4 of the MCTDHF non-linear derivative (CI
!>   with h2, orbitals with the projected correlation term) in imaginary time.
!> After each part the orbitals are Loewdin-orthonormalized with the
!> compensating CI transformation, so the many-body state is only changed by
!> the propagation itself (plain Gram-Schmidt would bias the fixed point), and
!> the CI vector is normalized.
!>
!> Initial orbitals either from orbsInit (program.lmaxIn < 0) or embedded from
!> orbXX_YY.in written on a grid with smaller lmax (program.lmaxIn >= 0, same
!> radial grid and grid.ylm.mmax); the CI vector from coeffsInit. The
!> imaginary time step is reduced (dtauFactor) whenever the energy change per
!> check drops below dtauSwitch, which removes the O(dtau^2) splitting bias
!> of the fixed point. The stiff non-linear part is integrated with
!> program.nSubStepsNonLin RK4 sub-steps per step, or, if
!> program.dtauNonLinMax is given, with as many as needed to keep every
!> sub-step at or below it (fewer as dtau shrinks).
!>
!> Output: orbXX_YY.in, coeffs.in (also written after every check, so a run
!>         can be restarted), rdm2BtXX_YY.in (TD-2RDM initial state),
!>         relax.dat (tau, dtau, energy)
program P_ImagTimeMctdhf
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Utils_RdmDiagonalize
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_IntegratorList_Sil

  implicit none

  real(R64)    :: tau, dtau, dtauMin, dtauFactor, dtauSwitch, convThresh, energy, energyOld, dtauNonLinMax
  integer(I32) :: lmaxIn, nStepsCheck, maxChecks, nSubStepsNonLin, iCheck, iStep, io
  real(R64), allocatable :: natOcc(:)
  complex(R64), allocatable :: rdm1(:, :), natOrbs(:, :)
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)

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
  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => ImagTimeDerivativeNonLin
  call IntegratorList_Fabricate(IntegratorListInput)

  call Say_Fabricate("program")
  lmaxIn = Json_Get("program.lmaxIn", -1)
  dtau = Json_Get("program.dtau", 0.05_R64)
  dtauMin = Json_Get("program.dtauMin", 0.005_R64)
  dtauFactor = Json_Get("program.dtauFactor", 0.5_R64)
  dtauSwitch = Json_Get("program.dtauSwitch", 1e-8_R64)
  convThresh = Json_Get("program.convThresh", 1e-11_R64)
  nStepsCheck = Json_Get("program.nStepsCheck", 20)
  maxChecks = Json_Get("program.maxChecks", 1000)
  nSubStepsNonLin = Json_Get("program.nSubStepsNonLin", 1)
  dtauNonLinMax = Json_Get("program.dtauNonLinMax", 0.0_R64)

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
  call IntegratorList_Setup

  select type (sil => integratorList(1) % e)
  type is (T_IntegratorList_E_Sil)
    if (.not. sil % imagTimeQ) error stop "integratorList(1) must have imagTimeQ = true"
  class default
    error stop "integratorList(1) must be a sil integrator"
  end select
  integratorList(1) % e % metricWeights = MetricWeights(Coeffs_nCoeffs)

  if (0 <= lmaxIn) call EmbedOrbitals(Orbs_orbs, lmaxIn)
  call LowdinKeepingState(Coeffs_coeffs, Orbs_orbs)
  call Coeffs_Normalize(Coeffs_coeffs)

  !==========================================
  call Say_Section("imaginary time propagation")
  !==========================================

  open (newunit=io, file="relax.dat", status="replace")

  energy = Method_GetEnergy(0.0_R64)
  write (*, '(A,F22.14)') " start energy: ", energy

  tau = 0.0_R64

  do iCheck = 1, maxChecks

    do iStep = 1, nStepsCheck
      call ImagTimeStep(dtau)
      tau = tau + dtau
    end do

    energyOld = energy
    energy = Method_GetEnergy(0.0_R64)

    write (io, '(2E14.6, F24.16)') tau, dtau, energy
    flush (io)

    ! checkpoint (restart with program.lmaxIn = grid lmax and coeffsInit.load)
    call Orbs_SaveOrbs(Orbs_orbs)
    call Coeffs_SaveCoeffs(Coeffs_coeffs)
    write (*, '(A,F10.3,A,ES9.2,A,F22.14,A,ES10.2)') " tau =", tau, "  dtau =", dtau, &
      "  energy =", energy, "  change =", energy - energyOld

    if (abs(energy - energyOld) < dtauSwitch .and. dtau > dtauMin * (1.0_R64 + 1e-12_R64)) then
      dtau = max(dtau * dtauFactor, dtauMin)
      cycle
    end if

    if (abs(energy - energyOld) < convThresh .and. dtau <= dtauMin * (1.0_R64 + 1e-12_R64)) exit

  end do

  close (io)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
  call RdmDiagonalize_Rdm1(natOcc, natOrbs, rdm1(1:Orbs_nOrbsInState, 1:Orbs_nOrbsInState))

  print *
  write (*, '(A,F22.14)') " final energy: ", energy
  write (*, '(A,*(ES12.4))') " natural occupations (spin up): ", natOcc

  !==========================================
  call Say_Section("save data")
  !==========================================

  call Orbs_SaveOrbs(Orbs_orbs)
  call Coeffs_SaveCoeffs(Coeffs_coeffs)
  call Coeffs_SaveTwoRdm(Coeffs_coeffs)

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> One Strang step exp(-h dtau/2) exp(-W dtau) exp(-h dtau/2) in imaginary
  !> time with state-preserving re-orthonormalization after each part.
  subroutine ImagTimeStep(dtauStep)
    real(R64), intent(in) :: dtauStep

    integer(I32) :: iSub, nSub

    call ImagTimeLin(0.5_R64 * dtauStep)
    call Reorthonormalize

    ! The non-linear part is stiff for weakly occupied orbitals (inverse
    ! 1-RDM), hence the explicit sub-steps
    nSub = nSubStepsNonLin
    if (dtauNonLinMax > 0.0_R64) nSub = max(1, ceiling(dtauStep / dtauNonLinMax * (1.0_R64 - 1e-12_R64)))
    do iSub = 1, nSub
      call integratorList(2) % e % Integrate(Method_state, 0.0_R64, dtauStep / nSub)
      call Reorthonormalize
    end do

    call ImagTimeLin(0.5_R64 * dtauStep)
    call Reorthonormalize

  end subroutine

  !> Loewdin orthonormalization of the orbitals with the compensating CI
  !> transformation, then normalization of the CI vector: the Loewdin step
  !> moves the orbital norms (exp(-eps dtau) per orbital) into the CI vector,
  !> and the non-linear derivative projects with 1 - |c><c|, which requires
  !> |c| = 1 (otherwise it amplifies, e.g. for Ne with eps(1s) = -50 a.u.).
  subroutine Reorthonormalize

    call LowdinKeepingState(Coeffs_coeffs, Orbs_orbs)
    call Coeffs_Normalize(Coeffs_coeffs)

  end subroutine

  !> Linear part exp(-h dtau) with integrator 1; with
  !> method.mb.orbBased.exactFieldFreeQ the field-free h0 is applied exactly
  !> around it (integrator 1 then only sees the vanishing gauge coupling).
  subroutine ImagTimeLin(dtauLin)
    real(R64), intent(in) :: dtauLin

    if (associated(Method_PropagateExactPart)) then
      call Method_PropagateExactPart(Method_state, 0.0_R64, 0.5_R64 * dtauLin, imagTimeQ_=.true.)
      call integratorList(1) % e % Integrate(Method_state, 0.0_R64, dtauLin)
      call Method_PropagateExactPart(Method_state, 0.0_R64, 0.5_R64 * dtauLin, imagTimeQ_=.true.)
    else
      call integratorList(1) % e % Integrate(Method_state, 0.0_R64, dtauLin)
    end if

  end subroutine

  !> Imaginary-time version of the non-linear split-step derivative:
  !> d/dtau = -i * (d/dt).
  subroutine ImagTimeDerivativeNonLin(dState, state, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in) :: time

    call Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin(dState, state, time)
    dState(:) = -IU * dState(:)

  end subroutine

  include '../common/BeCommon.inc'

end program
