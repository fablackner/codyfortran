!> ctest-requires: libxc
!> Real-time TDKS with a GGA (PBE) on a Ylm FEDVR grid: split step of the
!> exact field-free part and the Kohn-Sham potential at dt = 0.05 with the
!> spectral cut exactFieldFreeCut = 40 (Ec dt = 2 < pi).
!>
!> The splitting of h0 and a GGA potential is unstable once an h0 eigenvalue E
!> of the propagated space has E dt > pi (resonances of the exact rotation with
!> the GGA kernel, which grows like k^2): without the cut (FEDVR eigenvalues up
!> to ~1e5) the energy of this run grows to ~0.1 Ha, with Ec = 120 to ~1e-2.
!> The run starts from the PBE ground state projected below the cut (energy
!> shift ~1e-5 Ha) and checks energy conservation, the shift and that no
!> weight appears above the cut.
program T_He3d_04_TimePropagationTdksPbeCut
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Orbs
  use M_Method
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)
  real(R64) :: energyScf, energy0, energy, maxDrift, time, timeStep, endTime
  real(R64), allocatable :: eps(:)
  complex(R64), allocatable :: orbsCut(:, :)
  integer(I32) :: nScf, iStep, nSteps

  call Say_Hello
  call RunKsScfYlm(energyScf, eps, nScf, "test/simulationTests/He3d/T_He3d_04_TimePropagationTdksPbeCut.json", 40)

  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()
  call IntegratorList_Setup
  call Propagator_Setup

  ! the SCF orbitals have weight above the cut: start below it
  call Method_Mb_OrbBased_ApplySpectralCut(Orbs_orbs)
  call Orbs_Orthonormalize(Orbs_orbs)
  energy0 = Method_GetEnergy(0.0_R64)

  timeStep = Json_Get("program.timeStep", 0.05_R64)
  endTime = Json_Get("program.endTime", 8.0_R64)
  nSteps = nint(endTime / timeStep)
  time = 0.0_R64
  maxDrift = 0.0_R64
  do iStep = 1, nSteps
    call Propagator_Propagate(Method_state, time, time + timeStep)
    time = time + timeStep
    if (mod(iStep, 20) .eq. 0) then
      energy = Method_GetEnergy(time)
      maxDrift = max(maxDrift, abs(energy - energy0))
      write (*, '(1X, A, F8.3, A, ES11.3)') "t ", time, "  E - E0 ", energy - energy0
    end if
  end do

  orbsCut = Orbs_orbs
  call Method_Mb_OrbBased_ApplySpectralCut(orbsCut)

  print *
  print '(1X, A, F16.10, A, ES10.2, A, ES10.2, A, ES10.2)', "He PBE: E_SCF = ", energyScf, &
    ", cut shift ", energy0 - energyScf, ", max |E - E0| ", maxDrift, &
    ", weight above cut ", maxval(abs(orbsCut - Orbs_orbs))

  call check(error, maxDrift < 5e-6_R64)
  if (allocated(error)) error stop "T_He3d_04_TimePropagationTdksPbeCut failure: energy not conserved"
  call check(error, abs(energy0 - energyScf) < 5e-5_R64)
  if (allocated(error)) error stop "T_He3d_04_TimePropagationTdksPbeCut failure: energy shift of the cut"
  call check(error, maxval(abs(orbsCut - Orbs_orbs)) < 1e-10_R64)
  if (allocated(error)) error stop "T_He3d_04_TimePropagationTdksPbeCut failure: weight above the cut"

  call Say_Goodbye

contains

  include "../KsScfYlm.inc"

end program
