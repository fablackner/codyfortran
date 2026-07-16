!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Real-time MCTDHF propagation of Ne (3d, Ylm grid) in a strong z-polarized
!> femtosecond laser pulse (length gauge), computing the time-dependent
!> dipole moment.
!>
!> The correlated ground state is loaded from file (orb*.in, coeffs.in as
!> written by the mcscf ground-state app on the SAME grid — run the app in
!> the directory containing those files). The state is propagated with the
!> split-step (Strang) propagator: the stiff linear one-body part
!> (kinetic + centrifugal + Coulomb + laser) is integrated by the Short
!> Iterative Lanczos in the FEDVR-weighted metric, the nonlinear
!> mean-field/CI part by classical RK4.
!>
!> Output (appended per outputStep):
!>   dipole.dat : time, <z>
!>   laser.dat  : time, E(t)
!>   energy.dat : time, <H(t)>
program P_TimePropagationMctdhxLaser
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Actions
  use M_Utils_PrinterObservable
  use M_Utils_PrinterObservableYlm
  use M_Grid
  use M_Grid_Ylm
  use M_SysKinetic
  use M_SysPotential
  use M_SysPotential_Ylm_CoulombLaser
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator

  implicit none

  real(R64)    :: time, startTime, endTime, outputStep, timeStep
  real(R64)    :: startNorm, norm
  real(R64)    :: dipole, energy, field
  integer(I32) :: iStep, nTimeStepsPropagation, innerStep, iOrb
  integer(I32) :: io, istat
  real(R64), allocatable :: metricWeights(:)
  complex(R64), allocatable :: rdm1(:, :)
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
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Say_Fabricate("program")
  startTime = Json_Get("program.startTime", 0.0_R64)
  endTime = Json_Get("program.endTime", 1.0_R64)
  outputStep = Json_Get("program.outputStep", 1e-1_R64)
  nTimeStepsPropagation = Json_Get("program.nTimeStepsPropagation", 10)

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
  call Propagator_Setup

  ! The one-body operator is self-adjoint only in the FEDVR-weighted metric;
  ! supply the metric to the SIL integrator (integrator 1): unit weights for
  ! the CI coefficients, grid weights for each orbital.
  allocate (metricWeights(Coeffs_nCoeffs + Grid_nPoints * Orbs_nOrbsInState))
  metricWeights(1:Coeffs_nCoeffs) = 1.0_R64
  do iOrb = 1, Orbs_nOrbsInState
    metricWeights(Coeffs_nCoeffs + (iOrb - 1) * Grid_nPoints + 1: &
                  Coeffs_nCoeffs + iOrb * Grid_nPoints) = Grid_Ylm_weights(:)
  end do
  integratorList(1) % e % metricWeights = metricWeights

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  startNorm = MetricNorm(Method_state)
  energy = Method_GetEnergy(startTime)

  write (*, *)
  write (*, *) "start energy: ", energy
  write (*, *) "start norm:   ", startNorm
  write (*, *)

  time = startTime
  iStep = 0
  timeStep = outputStep / nTimeStepsPropagation

  do while (time < endTime - 0.5_R64 * timeStep)
    iStep = iStep + 1

    do innerStep = 1, nTimeStepsPropagation
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
    end do

    energy = Method_GetEnergy(time)
    field = SysPotential_Ylm_CoulombLaser_FieldAmplitude(time)
    norm = MetricNorm(Method_state)

    call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
    call PrinterObservableYlm_DumpDipole(rdm1, Orbs_orbs, "dipole.dat", .false., time, dipole_=dipole)
    call PrinterObservable_DumpEnergy("energy.dat", .false., time)

    open (newunit=io, file="laser.dat", status="unknown", position="append", iostat=istat)
    if (istat .eq. 0) then
      write (io, '(2E20.10E3)') time, field
      close (io)
    end if

    write (*, '(A,ES12.5E1,A,E20.10,A,E20.10,A,E14.6,A,F14.10)') &
      ' time = ', time, '    energy = ', energy, '    <z> = ', dipole, &
      '    E(t) = ', field, '    norm = ', norm
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "final energy: ", energy
  print *, "final dipole: ", dipole
  print *, "start norm:   ", startNorm
  print *, "final norm:   ", MetricNorm(Method_state)

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  function MetricNorm(state) result(res)
    complex(R64), intent(in), contiguous :: state(:)
    real(R64) :: res

    res = sqrt(sum(metricWeights(:) * abs(state(:))**2))

  end function

end program
