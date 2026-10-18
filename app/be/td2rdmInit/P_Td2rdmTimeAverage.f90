!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Initial state of the TD-2RDM propagation (Lackner et al., PRA 95, 033414
!> (2017), appendix C): the MCTDHF ground state is not stationary under the
!> TD-2RDM equations of motion once the 3-RDM is reconstructed, so the
!> residual excitations are projected out by a field-free TD-2RDM propagation
!> and the time average
!>   D12 = (1/T) int_0^T D12(t) dt .
!>
!> In the split-step gauge the orbitals rotate within their span (exp(-i h1 t))
!> and the 2-RDM compensates, so only the 2-RDM kernel
!>   D(x1 x2; x1' x2') = sum D(i1,i2,j1,j2) phi_i1(x1) phi_i2(x2) phi*_j1 phi*_j2
!> is (nearly) stationary. It is averaged in the fixed basis of the initial
!> orbitals phi0: with U = <phi0|phi(t)>, D_ref = (U x U) D (U x U)^+. The
!> small component of phi(t) outside span(phi0) is neglected. The averaged
!> 2-RDM is renormalized and stored together with phi0.
!>
!> Input:  orbXX_YY.in, rdm2BtXX_YY.in (relax/P_ImagTimeMctdhf)
!> Output: orbXX_YY.in, rdm2BtXX_YY.in (averaged),
!>         average.dat (time, energy, natural occupations along the propagation)
program P_Td2rdmTimeAverage
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_DataStorage
  use M_Utils_RdmDiagonalize
  use M_Grid
  use M_Grid_Ylm
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_TwoRdm
  use M_TwoRdmInit
  use M_Method
  use M_Method_Mb
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator
  use M_Reconstruction

  implicit none

  real(R64)    :: time, averagingTime, outputStep, timeStep, energy, weight
  integer(I32) :: nTimeSteps, innerStep, io
  complex(R64), allocatable :: orbs0(:, :), rdm2Average(:, :, :, :)
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
  call TwoRdm_Fabricate
  call TwoRdmInit_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()
  call Reconstruction_Fabricate

  call Say_Fabricate("program")
  averagingTime = Json_Get("program.averagingTime", 40.0_R64)
  outputStep = Json_Get("program.outputStep", 0.2_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 20)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call TwoRdm_Setup
  call TwoRdmInit_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup
  call Reconstruction_Setup

  ! The one-body operator is self-adjoint only in the FEDVR-weighted metric
  integratorList(1) % e % metricWeights = MetricWeights(TwoRdm_nElements)

  !==========================================
  call Say_Section("field-free propagation and time average")
  !==========================================

  open (newunit=io, file="average.dat", status="replace")

  time = 0.0_R64
  timeStep = outputStep / nTimeSteps

  orbs0 = Orbs_orbs

  ! trapezoidal rule over all time steps
  call AddToAverage(0.5_R64)
  call WriteProgress(time)

  do while (time < averagingTime - 0.5_R64 * timeStep)

    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep

      weight = 1.0_R64
      if (averagingTime - 0.5_R64 * timeStep <= time) weight = 0.5_R64
      call AddToAverage(weight)
    end do

    call WriteProgress(time)

  end do

  close (io)

  Orbs_orbs(:, :) = orbs0(:, :)
  call TwoRdm_FillTwoRdmStateFromRdm2(TwoRdm_twoRdmState, rdm2Average * (timeStep / time))
  call TwoRdm_Normalize(TwoRdm_twoRdmState)

  energy = Method_GetEnergy(0.0_R64)

  print *
  write (*, '(A,F22.14)') " energy of the averaged state: ", energy

  !==========================================
  call Say_Section("save data")
  !==========================================

  call Orbs_SaveOrbs(Orbs_orbs)
  call SaveTwoRdmBlocks

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Add weight * D(t), expressed in the initial orbital basis, to the average.
  subroutine AddToAverage(weight)
    real(R64), intent(in) :: weight

    complex(R64), allocatable :: rdm2(:, :, :, :), u(:, :), uFull(:, :)
    integer(I32) :: nO, nOS, i, j

    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    ! U(k,i) = <phi0_k|phi_i(t)>, same block for both spins when restricted
    allocate (u(nOS, nOS), uFull(nO, nO))
    do i = 1, nOS
      do j = 1, nOS
        u(j, i) = Grid_InnerProduct(orbs0(:, j), Orbs_orbs(:, i))
      end do
    end do
    uFull = 0.0_R64
    uFull(1:nOS, 1:nOS) = u
    if (Orbs_restrictedQ) uFull(nOS + 1:, nOS + 1:) = u

    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, TwoRdm_twoRdmState)
    call TransformRdm2(rdm2, uFull)

    if (.not. allocated(rdm2Average)) then
      allocate (rdm2Average, mold=rdm2)
      rdm2Average = 0.0_R64
    end if
    rdm2Average = rdm2Average + weight * rdm2

  end subroutine

  !> rdm2 <- sum U(k1,i1) U(k2,i2) rdm2(i1,i2,j1,j2) conjg(U(l1,j1)) conjg(U(l2,j2))
  subroutine TransformRdm2(rdm2, u)
    complex(R64), intent(inout), contiguous :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous    :: u(:, :)

    complex(R64), allocatable :: tmp(:, :, :, :)
    integer(I32) :: n, a, b, c

    n = size(u, 1)
    allocate (tmp, mold=rdm2)

    ! index 1 and 3 (creators/annihilators of particle 1), then 2 and 4
    do c = 1, n
      do b = 1, n
        do a = 1, n
          tmp(:, a, b, c) = matmul(u, rdm2(:, a, b, c))
        end do
      end do
    end do
    do c = 1, n
      do b = 1, n
        do a = 1, n
          rdm2(a, :, b, c) = matmul(u, tmp(a, :, b, c))
        end do
      end do
    end do
    do c = 1, n
      do b = 1, n
        do a = 1, n
          tmp(a, b, :, c) = matmul(conjg(u), rdm2(a, b, :, c))
        end do
      end do
    end do
    do c = 1, n
      do b = 1, n
        do a = 1, n
          rdm2(a, b, c, :) = matmul(conjg(u), tmp(a, b, c, :))
        end do
      end do
    end do

  end subroutine

  subroutine WriteProgress(t)
    real(R64), intent(in) :: t

    real(R64), allocatable :: natOcc(:)
    complex(R64), allocatable :: rdm1(:, :), natOrbs(:, :)

    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, TwoRdm_twoRdmState)
    call RdmDiagonalize_Rdm1(natOcc, natOrbs, rdm1(1:Orbs_nOrbsInState, 1:Orbs_nOrbsInState))

    energy = Method_GetEnergy(t)

    write (io, '(F10.4, F24.16, *(E16.6E3))') t, energy, natOcc
    flush (io)
    write (*, '(A,F9.3,A,F22.14)') " t =", t, "  energy =", energy

  end subroutine

  !> Write the body-type blocks rdm2BtXX_YY.in of the current 2-RDM (the
  !> format read by twoRdmInit.load and written by Coeffs_SaveTwoRdm).
  subroutine SaveTwoRdmBlocks
    complex(R64), allocatable :: rdm2(:, :, :, :), rdm2Bt(:, :, :, :)
    character(len=256) :: filename
    integer(I32) :: ibt1, ibt2

    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, TwoRdm_twoRdmState)

    do ibt1 = 1, Method_Mb_nBodyTypes
      do ibt2 = 1, Method_Mb_nBodyTypes
        write (filename, '("rdm2Bt",I2.2,"_",I2.2,".in")') ibt1, ibt2
        rdm2Bt = rdm2(Method_Mb_OrbBased_nOrbsStart(ibt1):Method_Mb_OrbBased_nOrbsEnd(ibt1), &
                      Method_Mb_OrbBased_nOrbsStart(ibt2):Method_Mb_OrbBased_nOrbsEnd(ibt2), &
                      Method_Mb_OrbBased_nOrbsStart(ibt1):Method_Mb_OrbBased_nOrbsEnd(ibt1), &
                      Method_Mb_OrbBased_nOrbsStart(ibt2):Method_Mb_OrbBased_nOrbsEnd(ibt2))
        call SaveData(rdm2Bt, trim(filename), storage_size(rdm2Bt))
      end do
    end do

  end subroutine

  include '../common/BeCommon.inc'

end program
