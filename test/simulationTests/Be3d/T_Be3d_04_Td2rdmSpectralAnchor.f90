!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Spectral anchor of the TD-2RDM reconstruction
!> (reconstruction.singlet.anchor.nTransitions): Be MCSCF state of T_Be3d_01
!> (orbitals, coefficients, 2-RDM, up-up-down 3-RDM), static orbitals
!> (td2rdmStaticOrbs, NY closure with ccQ).
!>
!> The test diagonalizes the CI Hamiltonian of the orbitals (Coeffs space,
!> total M = 0), takes the eigenstates k coupled to the ground state by the
!> one-body operator A = |2s><2p0| + |2p0><2s| (both spins), writes the
!> anchor files of c0 and these states (the 3-RDM of
!> c0 and the transition RDMs T2_k = <0|..|k>, T3_k), and checks:
!> - the closure at the reference, |D3_ref - R(D2_ref)| (the anchor hides
!>   closure errors, so this value pins the NY kernels);
!> - the reference is stationary (|dD2/dt| 8.0e-9 at the time of writing);
!> - the defining property: along D = D_ref + eps (t2 + t2^+), t2 = T2_k and
!>   t2 = i T2_k, the TD-2RDM time derivative equals the exact one,
!>   eps w_k (-i t2 + i t2^+), for every included excitation (relative error
!>   < 1e-5; 2e-8 at the time of writing);
!> - the last two again with the Valdemoro closure without ccQ
!>   (T_Be3d_04_Td2rdmSpectralAnchorValdemoro.json).
program T_Be3d_04_Td2rdmSpectralAnchor
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_DataStorage
  use M_Utils_LapackLib
  use M_Utils_Constants, only: CONE, IU
  use M_Utils_RdmTrace
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_TwoRdm
  use M_TwoRdm_Fermionic_Singlet
  use M_TwoRdmInit
  use M_Reconstruction
  use M_Reconstruction_Singlet_Anchor
  use M_Utils_RdmCumulants
  use M_Method_Mb
  use M_Method
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  complex(R64), allocatable :: c0(:), hCi(:, :), unit(:), evecs(:, :), x(:), a1(:, :), kets(:, :)
  complex(R64), allocatable :: trdm2(:, :, :, :, :), rdm2Ref(:, :, :, :)
  complex(R64), allocatable :: rdm1a(:, :), rdm3(:, :, :, :, :, :), rdm3Ref(:, :, :, :, :, :)
  complex(R64), allocatable :: dState(:)
  real(R64), allocatable :: evals(:), omega(:)
  integer(I32), allocatable :: selected(:)
  real(R64) :: amplitude, maxRelErr, maxRelErrNoCcQ, e0
  integer(I32) :: nC, nO, nOS, n, i, k, nT, nTransitions

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/simulationTests/Be3d/T_Be3d_04_Td2rdmSpectralAnchor.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call TwoRdm_Fabricate
  call TwoRdmInit_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  call Reconstruction_Fabricate

  call Say_Fabricate("program")
  amplitude = Json_Get("program.amplitude", 1e-4_R64)
  nTransitions = Json_Get("reconstruction.singlet.anchor.nTransitions", 0)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call TwoRdm_Setup
  call TwoRdmInit_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call Method_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup

  !==========================================
  call Say_Section("excited states and transition RDMs")
  !==========================================

  nC = Coeffs_nCoeffs
  nO = Method_Mb_OrbBased_nOrbsSum
  nOS = Orbs_nOrbsInState
  n = Method_Mb_OrbBased_nOrbs(1)

  allocate (c0(nC))
  call CoeffsInit_Initialize(c0)

  ! dense CI Hamiltonian of the static orbitals (h1, h2 filled by Method_Setup)
  allocate (hCi(nC, nC), unit(nC))
  do i = 1, nC
    unit = 0.0_R64
    unit(i) = 1.0_R64
    call Method_Mb_OrbBased_ApplyCiHamiltonian(hCi(:, i), unit, 0.0_R64)
  end do
  hCi = 0.5_R64 * (hCi + conjg(transpose(hCi)))
  call LapackLib_DiagonalizeGeneric(evals, evecs, hCi, .true.)
  e0 = real(dot_product(c0, matmul(hCi, c0)), R64)
  write (*, '(A,2F20.12)') " lowest CI eigenvalue, <c0|H|c0>: ", evals(1), e0

  ! A = |2s><2p0| + h.c. in both spin blocks; eigenstates coupled to c0
  allocate (a1(nO, nO), x(nC))
  a1 = 0.0_R64
  a1(2, 4) = 1.0_R64
  a1(4, 2) = 1.0_R64
  a1(nOS + 2, nOS + 4) = 1.0_R64
  a1(nOS + 4, nOS + 2) = 1.0_R64
  x = 0.0_R64
  call Coeffs_ApplyH1FillRdm1(c0, dCoeffs_=x, h1_=a1)
  allocate (selected(0))
  do k = 2, nC
    if (abs(dot_product(evecs(:, k), x))**2 > 1e-8_R64 .and. size(selected) < nTransitions) selected = [selected, k]
  end do
  nT = size(selected)
  if (nT < 1 .or. nT .ne. nTransitions) &
    error stop "T_Be3d_04: needs 1 <= reconstruction.singlet.anchor.nTransitions <= the eigenstates coupled to c0"
  omega = evals(selected) - e0
  write (*, '(A,*(F12.6))') " excitation energies: ", omega

  allocate (kets(nC, nT + 1), trdm2(n, n, n, n, nT))
  kets(:, 1) = c0
  kets(:, 2:) = evecs(:, selected)
  call Reconstruction_Singlet_Anchor_SaveReference(kets)
  call LoadData(trdm2, "trdm2Bt01_02.in", storage_size(trdm2))

  call Reconstruction_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2Ref, TwoRdm_twoRdmState)
  allocate (dState, mold=Method_state)

  ! the closure at the reference against the exact 3-RDM
  allocate (rdm3Ref(n, n, n, n, n, n))
  call LoadData(rdm3Ref, "rdm3Bt01_01_02.in", storage_size(rdm3Ref))
  call RdmTrace_FillRdm1aFromRdm2ab(rdm1a, rdm2Ref)
  call RdmCumulants_FillCumulFreePartM3aab(rdm3, rdm1a, rdm2Ref)
  call RdmCumulants_AddNyCumulantM3aab(rdm3, Method_Mb_nBodies(1), rdm1a, rdm2Ref)
  write (*, '(A,ES16.8)') " |D3_ref - R(D2_ref)|: ", norm2(abs(rdm3 - rdm3Ref))
  call check(error, norm2(abs(rdm3 - rdm3Ref)), 3.6209e-2_R64, rel=.true., thr=1e-4_R64)
  if (allocated(error)) &
    error stop "T_Be3d_04_Td2rdmSpectralAnchor failure: closure at the reference (update the value if T_Be3d_01 changed)"

  call CheckAnchor(maxRelErr)

  ! the same without ccQ (Valdemoro closure): the e_a enter unprojected
  call Json_LoadJsonFile(manualFileName_="test/simulationTests/Be3d/T_Be3d_04_Td2rdmSpectralAnchorValdemoro.json", &
                         relativeToProjectDirQ_=.true.)
  call Reconstruction_Fabricate
  call Reconstruction_Setup
  call CheckAnchor(maxRelErrNoCcQ)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  write (*, '(A,2ES12.4)') " max relative error (NY with ccQ, Valdemoro without): ", maxRelErr, maxRelErrNoCcQ

  call check(error, maxRelErr < 1e-5_R64 .and. maxRelErrNoCcQ < 1e-5_R64)
  if (allocated(error)) error stop "T_Be3d_04_Td2rdmSpectralAnchor failure: derivative along the transitions"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> The reference is stationary, and along v = t2 + t2^+ (t2 = T2_k, i T2_k)
  !> the linearized derivative (central difference) is the exact one; returns
  !> the largest relative error.
  subroutine CheckAnchor(maxRelErr)
    real(R64), intent(out) :: maxRelErr

    complex(R64), allocatable :: t2(:, :, :, :), t2Adj(:, :, :, :), dExact(:, :, :, :), dLin(:, :, :, :)
    real(R64) :: relErr
    integer(I32) :: k, b

    dLin = Derivative(rdm2Ref)
    write (*, '(A,ES12.4)') " |dD2/dt| at the reference: ", norm2(abs(dLin))
    call check(error, norm2(abs(dLin)) < 1e-6_R64)
    if (allocated(error)) error stop "T_Be3d_04_Td2rdmSpectralAnchor failure: reference not stationary"

    maxRelErr = 0.0_R64
    do k = 1, nT
      do b = 1, 2
        t2 = merge(CONE, IU, b .eq. 1) * trdm2(:, :, :, :, k)
        t2Adj = conjg(reshape(t2, [n, n, n, n], order=[3, 4, 1, 2]))
        dExact = amplitude * omega(k) * IU * (t2Adj - t2)
        dLin = 0.5_R64 * (Derivative(rdm2Ref + amplitude * (t2 + t2Adj)) - Derivative(rdm2Ref - amplitude * (t2 + t2Adj)))
        relErr = norm2(abs(dLin - dExact)) / norm2(abs(dExact))
        maxRelErr = max(maxRelErr, relErr)
        write (*, '(A,I3,A,I2,A,F12.6,A,ES12.4)') " transition", k, ", t2 = i^", b - 1, " T2:  w =", omega(k), &
          "  relative error of dD2/dt: ", relErr
      end do
    end do

  end subroutine

  !> dD2/dt of the TD-2RDM at the 2-RDM rdm2ab (static orbitals: the state is the 2-RDM)
  function Derivative(rdm2ab) result(dRdm2ab)
    complex(R64), intent(in) :: rdm2ab(:, :, :, :)
    complex(R64), allocatable :: dRdm2ab(:, :, :, :)

    call TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab(Method_state, rdm2ab)
    call Method_TimeDerivative(dState, Method_state, 0.0_R64)
    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(dRdm2ab, dState)

  end function

end program
