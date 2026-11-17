!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> coeffs.hubbard.noSpinSym against coeffs.generic: 3 + 2 fermions in 5
!> orbitals, random Hermitian h1 blocks (different for the two spins) and a
!> random on-site interaction sum_s U_s n_{s,up} n_{s,dn}. The two backends
!> number and phase the configurations differently, so the comparison is on
!> what does not depend on the basis: the spectrum of the dense Hamiltonian
!> (H1 + H2 applied to every unit vector) and the 1-RDM and 2-RDM of the
!> non-degenerate ground state. For coeffs.hubbard the test also checks that
!> H1 and H2 add to dCoeffs_, that the RDMs contract to the ground-state
!> energy (RdmObservables_Energy) and that the excitation-based 2-RDM blocks
!> (Coeffs_FillTransitionRdm2Bt through Coeffs_ApplyExcitation) equal rdm2.
program T_Coeffs_Hubbard_VsGeneric
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_LapackLib
  use M_Utils_RdmObservables
  use M_Grid
  use M_Method
  use M_Method_Mb
  use M_Method_Mb_OrbBased
  use M_ConfigList
  use M_Coeffs
  use testdrive, only: check, error_type

  implicit none

  character(len=*), parameter :: dir = "test/componentTests/Coeffs/Hubbard/"
  type(error_type), allocatable :: error
  complex(R64), allocatable :: h1(:, :), h2(:, :, :, :), hUp(:, :), hDn(:, :), hMat(:, :), evecs(:, :)
  complex(R64), allocatable :: rdm1Generic(:, :), rdm2Generic(:, :, :, :), rdm1Hubbard(:, :), rdm2Hubbard(:, :, :, :)
  complex(R64), allocatable :: ground(:), y(:), d(:), t2(:, :, :, :, :)
  real(R64), allocatable :: evalsGeneric(:), evalsHubbard(:), u(:), yRe(:)
  real(R64) :: errSpectrum, errRdm1, errRdm2, errEnergy, errAccumulate, errBlocks, errHermitian, gap
  integer(I32) :: n, b1, b2

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=dir//"T_Coeffs_Hubbard_VsGeneric.json", relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Method_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate

  call Grid_Setup
  call ConfigList_Setup
  call Coeffs_Setup

  !==========================================
  call Say_Section("random h1 of both spins and on-site interaction")
  !==========================================

  n = Method_Mb_OrbBased_nOrbs(1)
  call FillRandomHermitian(hUp, n)
  call FillRandomHermitian(hDn, n)
  allocate (u(n))
  call random_number(u)
  u = 0.5_R64 + 2.0_R64 * u
  call FillLatticeHamiltonian(h1, h2, hUp, hDn, u)

  !==========================================
  call Say_Section("coeffs.generic")
  !==========================================

  call FillDenseHamiltonian(hMat, h1, h2)
  call LapackLib_DiagonalizeGeneric(evalsGeneric, evecs, hMat, .true.)
  ground = evecs(:, 1)
  call Coeffs_ApplyH1FillRdm1(ground, rdm1_=rdm1Generic)
  call Coeffs_ApplyH2FillRdm2(ground, rdm2_=rdm2Generic)

  !==========================================
  call Say_Section("coeffs.hubbard.noSpinSym")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=dir//"T_Coeffs_Hubbard_VsGenericHubbard.json", relativeToProjectDirQ_=.true.)
  call Coeffs_Fabricate
  call Coeffs_Setup

  call FillDenseHamiltonian(hMat, h1, h2)
  errHermitian = maxval(abs(hMat - conjg(transpose(hMat))))
  call LapackLib_DiagonalizeGeneric(evalsHubbard, evecs, hMat, .true.)
  ground = evecs(:, 1)
  gap = evalsHubbard(2) - evalsHubbard(1)
  call Coeffs_ApplyH1FillRdm1(ground, rdm1_=rdm1Hubbard)
  call Coeffs_ApplyH2FillRdm2(ground, rdm2_=rdm2Hubbard)
  errEnergy = abs(RdmObservables_Energy(rdm1Hubbard, h1, rdm2Hubbard, h2) - evalsHubbard(1))

  ! H1 and H2 add to what dCoeffs_ holds
  allocate (yRe(Coeffs_nCoeffs))
  call random_number(yRe)
  y = cmplx(yRe, 1.0_R64 - yRe, kind=R64)
  d = y
  call Coeffs_ApplyH1FillRdm1(ground, dCoeffs_=d, h1_=h1)
  call Coeffs_ApplyH2FillRdm2(ground, dCoeffs_=d, h2_=h2)
  errAccumulate = maxval(abs(d - y - evalsHubbard(1) * ground))

  ! 2-RDM blocks of the excitation route (Coeffs_ApplyExcitation)
  errBlocks = 0.0_R64
  do b2 = 1, 2
    do b1 = 1, 2
      call Coeffs_FillTransitionRdm2Bt(t2, ground, reshape(ground, [Coeffs_nCoeffs, 1]), b1, b2)
      associate (s1 => Method_Mb_OrbBased_nOrbsStart(b1), e1 => Method_Mb_OrbBased_nOrbsEnd(b1), &
                 s2 => Method_Mb_OrbBased_nOrbsStart(b2), e2 => Method_Mb_OrbBased_nOrbsEnd(b2))
        errBlocks = max(errBlocks, maxval(abs(t2(:, :, :, :, 1) - rdm2Hubbard(s1:e1, s2:e2, s1:e1, s2:e2))))
      end associate
    end do
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  errSpectrum = maxval(abs(evalsHubbard - evalsGeneric))
  errRdm1 = maxval(abs(rdm1Hubbard - rdm1Generic))
  errRdm2 = maxval(abs(rdm2Hubbard - rdm2Generic))

  print *
  print *, "coefficients:                    ", Coeffs_nCoeffs
  print *, "ground-state energy, gap:        ", evalsHubbard(1), gap
  print *, "max |H - H^+| (hubbard):         ", errHermitian
  print *, "max |spectrum hubbard - generic|:", errSpectrum
  print *, "max |rdm1 hubbard - generic|:    ", errRdm1
  print *, "max |rdm2 hubbard - generic|:    ", errRdm2
  print *, "|E(RDMs) - E0|:                  ", errEnergy
  print *, "max |d - y - H c| (accumulate):  ", errAccumulate
  print *, "max |2-RDM blocks - rdm2|:       ", errBlocks

  call check(error, errHermitian < 1e-13_R64)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_VsGeneric failure: H not Hermitian"
  call check(error, errSpectrum < 1e-11_R64)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_VsGeneric failure: spectrum differs from coeffs.generic"
  ! the RDM comparison needs a non-degenerate ground state
  call check(error, gap > 1e-3_R64)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_VsGeneric failure: degenerate ground state"
  call check(error, errRdm1 < 1e-11_R64 .and. errRdm2 < 1e-11_R64)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_VsGeneric failure: RDMs differ from coeffs.generic"
  call check(error, errEnergy < 1e-11_R64)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_VsGeneric failure: energy of the RDMs"
  call check(error, errAccumulate < 1e-11_R64)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_VsGeneric failure: H1, H2 do not add to dCoeffs_"
  call check(error, errBlocks < 1e-12_R64)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_VsGeneric failure: excitation-based 2-RDM blocks"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'CoeffsHubbardDense.inc'

end program
