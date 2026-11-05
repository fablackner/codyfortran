!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ne strong-field MCTDHF (9 restricted orbitals 1s 2s 2p 3s 3p, full CI in
!> the total-M = 0 sector) in a z-polarized 800 nm two-cycle velocity-gauge
!> pulse of 5e14 W/cm^2, started from the ground state of groundstate/mctdhf
!> (read from orbsInit.load.dir and coeffsInit.load.dir). It is relaxed at
!> lmax = 1, which is exact for the 1S state (its orbitals keep a pure l), and
!> embedded into the lmax = 40 propagation grid (orbsInit.load.lmaxIn). The
!> pulse lasts t = 0 .. 220.6 a.u. (88 240 steps,
!> ~1.6-2.4 s each on 4 threads).
!>
!> Split-step (Strang) propagation: the field-free h0 exactly (spectral cut
!> 1000 a.u.) around a Short Iterative Lanczos step for the A.p coupling in
!> the FEDVR metric; the interaction part by RK4. A radial cos^(1/order) mask
!> absorbs outgoing flux after every time step; with absorber.keepOrthonormalQ
!> the orbitals are then re-orthonormalized and the absorbed norm moves into
!> the CI vector (Method_Mb_OrbBased_OrthonormalizeKeepingState).
!>
!> Checkpoints: with program.saveStep > 0 (a multiple of program.outputStep)
!> the orbitals and the CI vector are saved (orbXX_YY.in, coeffs.in) after
!> every saveStep. A run is continued from the last checkpoint in the same
!> directory with orbsInit.load and coeffsInit.load without dir and lmaxIn
!> and program.startTime set to the checkpoint time;
!> observables.dat and natocc.dat then keep their lines before startTime and
!> are continued.
!>
!> Output (one line per program.outputStep):
!>   observables.dat : time, E(t), A(t), <z>, dipole acceleration
!>                     -Z<z/r^3>, electrons left in the box, energy
!>   natocc.dat      : time, natural occupation numbers (spin up)
!> and at the end the HHG spectrum of the whole observables.dat:
!>   hhg.dat         : harmonic order of sysGauge.ylm.velocity.omega,
!>                     acceleration form |a(w)|^2, length form w^4 |z(w)|^2
!>                     (Hann window, M_Utils_HhgSpectrum, up to
!>                     program.hhgMaxOrder, default 300)
program P_LaserMctdhfVelocity
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_HhgSpectrum
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysGauge
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
  use M_Absorber

  implicit none

  real(R64)    :: time, startTime, endTime, outputStep, timeStep, charge, saveStep
  integer(I32) :: nTimeSteps, innerStep, nOutputs, nOutputsPerSave
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile()

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysGauge_Fabricate
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
  call Absorber_Fabricate

  call Say_Fabricate("program")
  startTime = Json_Get("program.startTime", 0.0_R64)
  endTime = Json_Get("program.endTime", 1.0_R64)
  outputStep = Json_Get("program.outputStep", 0.2_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 20)
  saveStep = Json_Get("program.saveStep", 0.0_R64)
  charge = Json_Get("sysPotential.ylm.coulomb.charge", 4.0_R64)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysGauge_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup
  call Absorber_Setup

  ! The one-body operator is self-adjoint only in the FEDVR-weighted metric
  integratorList(1) % e % metricWeights = MetricWeights(Coeffs_nCoeffs)

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  time = startTime
  timeStep = outputStep / nTimeSteps
  nOutputsPerSave = nint(saveStep / outputStep)
  if (saveStep > 0.0_R64 .and. abs(nOutputsPerSave * outputStep - saveStep) > 1e-9_R64 * saveStep) &
    error stop "program.saveStep must be a multiple of program.outputStep"

  if (startTime > 0.0_R64) call TruncateLaserObservables(startTime - 0.5_R64 * timeStep)
  call WriteLaserObservables(time, charge, Method_GetEnergy(time), startTime <= 0.0_R64)

  nOutputs = 0
  do while (time < endTime - 0.5_R64 * timeStep)

    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
      call Absorber_ApplyAbsorber(Orbs_orbs)
      if (Absorber_keepOrthonormalQ) call Method_Mb_OrbBased_OrthonormalizeKeepingState(Method_state)
    end do

    call WriteLaserObservables(time, charge, Method_GetEnergy(time), .false.)

    nOutputs = nOutputs + 1
    if (nOutputsPerSave > 0) then
      if (mod(nOutputs, nOutputsPerSave) .eq. 0) then
        call Orbs_SaveOrbs(Orbs_orbs)
        call Coeffs_SaveCoeffs(Coeffs_coeffs)
        write (*, '(A,F9.3)') " checkpoint saved at t =", time
      end if
    end if

  end do

  !==========================================
  call Say_Section("save data")
  !==========================================

  call Orbs_SaveOrbs(Orbs_orbs)
  call Coeffs_SaveCoeffs(Coeffs_coeffs)

  !==========================================
  call Say_Section("hhg spectrum")
  !==========================================

  call HhgSpectrum_WriteFromTable("hhg.dat", "observables.dat", Json_Get("sysGauge.ylm.velocity.omega", 0.057_R64), &
                                  1, 4, colAcceleration_=5, maxOrder_=Json_Get("program.hhgMaxOrder", 300.0_R64))

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> FEDVR metric weights of a packed state [many-body part (nHead entries,
  !> unit weights) | orbitals (Grid_Ylm_weights each)] for the SIL integrator.
  function MetricWeights(nHead) result(res)
    use M_Grid
    use M_Grid_Ylm
    use M_Orbs

    integer(I32), intent(in) :: nHead
    real(R64), allocatable :: res(:)

    integer(I32) :: iOrb

    allocate (res(nHead + Grid_nPoints * Orbs_nOrbsInState))
    res(1:nHead) = 1.0_R64
    do iOrb = 1, Orbs_nOrbsInState
      res(nHead + (iOrb - 1) * Grid_nPoints + 1:nHead + iOrb * Grid_nPoints) = Grid_Ylm_weights(:)
    end do

  end function

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> The Hermitian matrix U diag(f) U^H of eigenvectors U and values f.
  function HermitianFunction(evecs, f) result(res)
    complex(R64), intent(in), contiguous :: evecs(:, :)
    real(R64), intent(in), contiguous    :: f(:)
    complex(R64), allocatable :: res(:, :)

    res = matmul(evecs, spread(f, 2, size(f)) * conjg(transpose(evecs)))

  end function

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Append one line to observables.dat and natocc.dat (replaced if newFilesQ):
  !>   observables.dat: time, E(t), A(t), <z>, -Z<z/r^3>, electrons in the box,
  !>                    energy
  !>   natocc.dat:      time, natural occupations of the spin-up 1-RDM
  !> Absorbed orbitals are no longer orthonormal: the electron number is
  !> Tr(rho S) and the occupations are the eigenvalues of S^1/2 rho S^1/2 with
  !> the orbital overlap matrix S (i.e. of the density matrix inside the box).
  !> The 1-RDM is taken from Method_state via the method-agnostic
  !> Method_Mb_OrbBased_FillRdm1 (MCTDHF and TD-2RDM alike).
  subroutine WriteLaserObservables(time, charge, energy, newFilesQ)
    use M_Utils_LapackLib
    use M_Utils_PrinterObservableYlm
    use M_Grid
    use M_Orbs
    use M_SysGauge_Ylm_Velocity
    use M_Method
    use M_Method_Mb_OrbBased

    real(R64), intent(in) :: time
    real(R64), intent(in) :: charge
    real(R64), intent(in) :: energy
    logical, intent(in) :: newFilesQ

    real(R64) :: dipole, acceleration, nElectrons
    real(R64), allocatable :: evals(:), occ(:)
    complex(R64), allocatable :: rdm1(:, :), s(:, :), evecs(:, :), sSqrt(:, :), rhoUp(:, :)
    integer(I32) :: nOS, i, j, io
    character(len=7) :: openStatus

    nOS = Orbs_nOrbsInState

    call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)

    call PrinterObservableYlm_DumpDipole(rdm1, Orbs_orbs, "", .false., time, dipole_=dipole)
    call PrinterObservableYlm_DumpDipoleAcceleration(rdm1, Orbs_orbs, charge, "", .false., time, &
                                                     acceleration_=acceleration)

    ! Overlap of the stored orbitals and spin-up block of the 1-RDM
    allocate (s(nOS, nOS))
    do j = 1, nOS
      do i = 1, nOS
        s(i, j) = Grid_InnerProduct(Orbs_orbs(:, i), Orbs_orbs(:, j))
      end do
    end do
    rhoUp = rdm1(1:nOS, 1:nOS)

    ! rdm1(i,j) = <a+_j a_i>: N = sum_ij rdm1(i,j) <phi_j|phi_i> per stored set
    nElectrons = real(sum(rhoUp * transpose(s)), kind=R64) * (size(rdm1, 1) / nOS)

    call LapackLib_DiagonalizeGeneric(evals, evecs, s, .true.)
    sSqrt = HermitianFunction(evecs, sqrt(max(evals, 0.0_R64)))
    call LapackLib_DiagonalizeGeneric(occ, evecs, matmul(conjg(transpose(sSqrt)), matmul(rhoUp, sSqrt)), .false.)

    openStatus = "unknown"
    if (newFilesQ) openStatus = "replace"

    open (newunit=io, file="observables.dat", status=trim(openStatus), position="append")
    if (newFilesQ) write (io, '(A)') "# time  E(t)  A(t)  <z>  -Z<z/r^3>  N_box  energy"
    write (io, '(7E24.14E3)') time, SysGauge_Ylm_Velocity_FieldAmplitude(time), &
      SysGauge_Ylm_Velocity_VectorPotentialAmplitude(time), dipole, acceleration, nElectrons, energy
    close (io)

    open (newunit=io, file="natocc.dat", status=trim(openStatus), position="append")
    write (io, '(*(E16.6E3))') time, occ(size(occ):1:-1)
    close (io)

    write (*, '(A,F9.3,A,ES12.4,A,ES14.6,A,ES14.6,A,F9.6,A,F18.10)') " t =", time, "  E =", &
      SysGauge_Ylm_Velocity_FieldAmplitude(time), "  <z> =", dipole, "  a =", acceleration, &
      "  N =", nElectrons, "  energy =", energy

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Remove the lines with time > tCut from observables.dat and natocc.dat (if
  !> present) before a run is continued from a checkpoint; comment lines stay.
  subroutine TruncateLaserObservables(tCut)
    real(R64), intent(in) :: tCut

    call TruncateTimeSeries("observables.dat", tCut)
    call TruncateTimeSeries("natocc.dat", tCut)

  end subroutine

  !> Keep the lines of a time series file up to time tCut (first column).
  subroutine TruncateTimeSeries(filename, tCut)
    character(len=*), intent(in) :: filename
    real(R64), intent(in) :: tCut

    character(len=4096), allocatable :: lines(:)
    character(len=4096) :: line
    real(R64) :: t
    integer(I32) :: io, ios, i
    logical :: existsQ

    inquire (file=filename, exist=existsQ)
    if (.not. existsQ) return

    allocate (lines(0))
    open (newunit=io, file=filename, status="old", action="read")
    do
      read (io, '(A)', iostat=ios) line
      if (ios .ne. 0) exit
      if (line(1:1) .ne. "#") then
        read (line, *) t
        if (t > tCut) exit
      end if
      lines = [lines, line]
    end do
    close (io)

    open (newunit=io, file=filename, status="replace", action="write")
    do i = 1, size(lines)
      write (io, '(A)') trim(lines(i))
    end do
    close (io)

  end subroutine

end program
