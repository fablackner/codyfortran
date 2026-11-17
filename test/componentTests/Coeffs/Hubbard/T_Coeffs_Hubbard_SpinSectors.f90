!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The spin-exchange sectors of coeffs.hubbard against the product space:
!> 2 + 2 fermions in 4 orbitals, a random Hermitian h1 that is the same for
!> both spins and a random on-site interaction. The spectra of the dense
!> Hamiltonians of plusSpinSym (21 coefficients) and minusSpinSym (15) make
!> up the spectrum of noSpinSym (36), and the lowest state of each sector has
!> the 1-RDM and 2-RDM of the noSpinSym eigenstate of the same energy (the
!> sector vectors are expanded to the product space for the RDMs).
program T_Coeffs_Hubbard_SpinSectors
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_LapackLib
  use M_Grid
  use M_Method
  use M_Method_Mb_OrbBased
  use M_Coeffs
  use testdrive, only: check, error_type

  implicit none

  character(len=*), parameter :: dir = "test/componentTests/Coeffs/Hubbard/"
  type(error_type), allocatable :: error
  complex(R64), allocatable :: h1(:, :), h2(:, :, :, :), hs(:, :), hMat(:, :), evecs(:, :)
  complex(R64), allocatable :: rdm1Plus(:, :), rdm2Plus(:, :, :, :), rdm1Minus(:, :), rdm2Minus(:, :, :, :)
  complex(R64), allocatable :: rdm1Full(:, :), rdm2Full(:, :, :, :)
  real(R64), allocatable :: evalsPlus(:), evalsMinus(:), evalsFull(:), evalsUnion(:), u(:)
  real(R64) :: errSpectrum, errHermitian, errRdmPlus, errRdmMinus, errTrace
  integer(I32) :: n, nPlus, nMinus, i

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=dir//"T_Coeffs_Hubbard_SpinSectorsPlus.json", relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Method_Fabricate
  call Coeffs_Fabricate

  call Grid_Setup
  call Coeffs_Setup

  !==========================================
  call Say_Section("random spin-independent h1 and on-site interaction")
  !==========================================

  n = Method_Mb_OrbBased_nOrbs(1)
  call FillRandomHermitian(hs, n)
  allocate (u(n))
  call random_number(u)
  u = 0.5_R64 + 2.0_R64 * u
  call FillLatticeHamiltonian(h1, h2, hs, hs, u)

  !==========================================
  call Say_Section("plusSpinSym")
  !==========================================

  nPlus = Coeffs_nCoeffs
  call FillDenseHamiltonian(hMat, h1, h2)
  errHermitian = maxval(abs(hMat - conjg(transpose(hMat))))
  call LapackLib_DiagonalizeGeneric(evalsPlus, evecs, hMat, .true.)
  call Coeffs_ApplyH1FillRdm1(evecs(:, 1), rdm1_=rdm1Plus)
  call Coeffs_ApplyH2FillRdm2(evecs(:, 1), rdm2_=rdm2Plus)
  errTrace = abs(sum([(rdm1Plus(i, i), i=1, size(rdm1Plus, 1))]) - 4.0_R64)

  !==========================================
  call Say_Section("minusSpinSym")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=dir//"T_Coeffs_Hubbard_SpinSectorsMinus.json", relativeToProjectDirQ_=.true.)
  call Coeffs_Fabricate
  call Coeffs_Setup

  nMinus = Coeffs_nCoeffs
  call FillDenseHamiltonian(hMat, h1, h2)
  errHermitian = max(errHermitian, maxval(abs(hMat - conjg(transpose(hMat)))))
  call LapackLib_DiagonalizeGeneric(evalsMinus, evecs, hMat, .true.)
  call Coeffs_ApplyH1FillRdm1(evecs(:, 1), rdm1_=rdm1Minus)
  call Coeffs_ApplyH2FillRdm2(evecs(:, 1), rdm2_=rdm2Minus)

  !==========================================
  call Say_Section("noSpinSym")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=dir//"T_Coeffs_Hubbard_SpinSectors.json", relativeToProjectDirQ_=.true.)
  call Coeffs_Fabricate
  call Coeffs_Setup

  call FillDenseHamiltonian(hMat, h1, h2)
  call LapackLib_DiagonalizeGeneric(evalsFull, evecs, hMat, .true.)

  evalsUnion = [evalsPlus, evalsMinus]
  call Sort(evalsUnion)
  errSpectrum = maxval(abs(evalsUnion - evalsFull))

  call FillRdmsOfEnergy(rdm1Full, rdm2Full, evalsPlus(1))
  errRdmPlus = max(maxval(abs(rdm1Plus - rdm1Full)), maxval(abs(rdm2Plus - rdm2Full)))
  call FillRdmsOfEnergy(rdm1Full, rdm2Full, evalsMinus(1))
  errRdmMinus = max(maxval(abs(rdm1Minus - rdm1Full)), maxval(abs(rdm2Minus - rdm2Full)))

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "coefficients plus, minus, product:  ", nPlus, nMinus, Coeffs_nCoeffs
  print *, "lowest energies plus, minus:        ", evalsPlus(1), evalsMinus(1)
  print *, "max |H - H^+| of the sectors:       ", errHermitian
  print *, "max |spectra plus + minus - full|:  ", errSpectrum
  print *, "|tr rdm1 - N| (plus):               ", errTrace
  print *, "max |RDMs plus - product space|:    ", errRdmPlus
  print *, "max |RDMs minus - product space|:   ", errRdmMinus

  call check(error, nPlus + nMinus, Coeffs_nCoeffs)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_SpinSectors failure: sector dimensions"
  call check(error, errHermitian < 1e-13_R64)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_SpinSectors failure: sector Hamiltonian not Hermitian"
  call check(error, errSpectrum < 1e-11_R64)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_SpinSectors failure: sector spectra differ from the product space"
  call check(error, errTrace < 1e-12_R64)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_SpinSectors failure: particle number of a sector state"
  call check(error, errRdmPlus < 1e-11_R64 .and. errRdmMinus < 1e-11_R64)
  if (allocated(error)) error stop "T_Coeffs_Hubbard_SpinSectors failure: RDMs of a sector state"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> RDMs of the noSpinSym eigenvector (evecs, evalsFull) of the energy e,
  !> which must be non-degenerate.
  subroutine FillRdmsOfEnergy(rdm1, rdm2, e)
    complex(R64), intent(out), allocatable :: rdm1(:, :)
    complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
    real(R64), intent(in) :: e

    integer(I32) :: j

    j = minloc(abs(evalsFull - e), 1)
    if (count(abs(evalsFull - e) < 1e-6_R64) .ne. 1) error stop "T_Coeffs_Hubbard_SpinSectors: degenerate eigenvalue"
    call Coeffs_ApplyH1FillRdm1(evecs(:, j), rdm1_=rdm1)
    call Coeffs_ApplyH2FillRdm2(evecs(:, j), rdm2_=rdm2)

  end subroutine

  !> Ascending insertion sort.
  subroutine Sort(a)
    real(R64), intent(inout), contiguous :: a(:)

    integer(I32) :: i, j
    real(R64) :: v

    do i = 2, size(a)
      v = a(i)
      j = i - 1
      do while (j >= 1)
        if (a(j) <= v) exit
        a(j + 1) = a(j)
        j = j - 1
      end do
      a(j + 1) = v
    end do

  end subroutine

  include 'CoeffsHubbardDense.inc'

end program
