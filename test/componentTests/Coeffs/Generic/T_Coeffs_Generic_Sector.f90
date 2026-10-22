program T_Coeffs_Generic_Sector
  !-----------------------------------------------------------------------------
  ! CI restricted to the declared total-M sector (coeffs.generic.totalM = 1)
  ! against the full space: three body types (2 and 3 fermions in 5 orbitals
  ! of m = 0 0 -1 0 1, 2 bosons in 3 orbitals of m = -1 0 1), random h1, h2
  ! (also with Delta M /= 0 entries) and a random vector of total M = 1.
  ! Coeffs_ApplyH1FillRdm1 and Coeffs_ApplyH2FillRdm2 of the sector give the
  ! full rdm1, rdm2 and the in-sector part of the full H|c>, and leave the
  ! rows outside the sector exactly zero. Coeffs_ProjectOnSector removes a
  ! small admixture of other M exactly. The M of a configuration is taken
  ! here from the number operators a+_i a_i (ConfigList ExciteConfiguration).
  !-----------------------------------------------------------------------------
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Method
  use M_Method_Mb
  use M_Method_Mb_OrbBased
  use M_ConfigList
  use M_Coeffs
  use testdrive, only: check, error_type

  implicit none

  character(len=*), parameter :: dir = "test/componentTests/Coeffs/Generic/"

  type(error_type), allocatable :: error
  complex(R64), allocatable :: coeffs(:), h1(:, :), h2(:, :, :, :)
  complex(R64), allocatable :: d1Sector(:), d2Sector(:), d1Full(:), d2Full(:)
  complex(R64), allocatable :: rdm1Sector(:, :), rdm2Sector(:, :, :, :), rdm1Full(:, :), rdm2Full(:, :, :, :)
  complex(R64), allocatable :: projected(:)
  logical, allocatable :: inSectorQ(:)
  integer(I32) :: nO, nNonZeroOutside, nChangedProjected
  real(R64) :: errRdm1, errRdm2, errH1, errH2

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=dir//"T_Coeffs_Generic_Sector.json", relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Method_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate

  call Grid_Setup
  call ConfigList_Setup
  call Coeffs_Setup

  !==========================================
  call Say_Section("random vector of the sector, h1, h2")
  !==========================================

  nO = Method_Mb_OrbBased_nOrbsSum
  inSectorQ = InSector(Json_Get("coeffs.generic.totalM", 0))

  h1 = reshape(RandomVector(nO**2), [nO, nO])
  h2 = reshape(RandomVector(nO**4), [nO, nO, nO, nO])
  coeffs = RandomVector(Coeffs_nCoeffs)
  where (.not. inSectorQ) coeffs = 0.0_R64
  call Coeffs_Normalize(coeffs)

  !==========================================
  call Say_Section("sector")
  !==========================================

  call Apply(d1Sector, d2Sector, rdm1Sector, rdm2Sector)

  ! Admixture of other M of relative norm ~1e-12, below the error threshold
  projected = coeffs + merge((0.0_R64, 0.0_R64), 1e-6_R64 * RandomVector(Coeffs_nCoeffs), inSectorQ)
  call Coeffs_ProjectOnSector(projected)
  nChangedProjected = count(projected .ne. merge(coeffs, (0.0_R64, 0.0_R64), inSectorQ))

  !==========================================
  call Say_Section("full space")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=dir//"T_Coeffs_Generic_SectorFull.json", relativeToProjectDirQ_=.true.)
  call Coeffs_Fabricate
  call Coeffs_Setup

  call Apply(d1Full, d2Full, rdm1Full, rdm2Full)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  errRdm1 = maxval(abs(rdm1Sector - rdm1Full))
  errRdm2 = maxval(abs(rdm2Sector - rdm2Full))
  errH1 = maxval(abs(d1Sector - d1Full), mask=inSectorQ)
  errH2 = maxval(abs(d2Sector - d2Full), mask=inSectorQ)
  nNonZeroOutside = count(.not. inSectorQ .and. (d1Sector .ne. 0.0_R64 .or. d2Sector .ne. 0.0_R64))

  print *
  print *, "sector configurations:          ", count(inSectorQ), " of ", Coeffs_nCoeffs
  print *, "max |rdm1 sector - full|:       ", errRdm1
  print *, "max |rdm2 sector - full|:       ", errRdm2
  print *, "max |H1 c sector - P H1 c full|:", errH1
  print *, "max |H2 c sector - P H2 c full|:", errH2
  print *, "max |P_out H2 c full|:          ", maxval(abs(d2Full), mask=.not. inSectorQ)
  print *, "nonzero H c outside the sector: ", nNonZeroOutside
  print *, "entries changed by projection:  ", nChangedProjected

  call check(error, errRdm1 < 1e-14_R64 .and. errRdm2 < 1e-14_R64)
  if (allocated(error)) error stop "T_Coeffs_Generic_Sector failure: RDMs of sector and full space differ"

  call check(error, errH1 < 1e-14_R64 .and. errH2 < 1e-14_R64)
  if (allocated(error)) error stop "T_Coeffs_Generic_Sector failure: H|c> of sector and full space differ"

  ! Rows outside the sector are not touched: exact zeros
  call check(error, nNonZeroOutside, 0)
  if (allocated(error)) error stop "T_Coeffs_Generic_Sector failure: H|c> outside the sector"

  ! Exactly the sector part of the vector survives the projection
  call check(error, nChangedProjected, 0)
  if (allocated(error)) error stop "T_Coeffs_Generic_Sector failure: Coeffs_ProjectOnSector"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> H1|c>, H2|c> (accumulated on zero), rdm1 and rdm2 of the bound backend.
  subroutine Apply(d1, d2, rdm1, rdm2)
    complex(R64), allocatable, intent(out) :: d1(:)
    complex(R64), allocatable, intent(out) :: d2(:)
    complex(R64), allocatable, intent(out) :: rdm1(:, :)
    complex(R64), allocatable, intent(out) :: rdm2(:, :, :, :)

    allocate (d1(Coeffs_nCoeffs), d2(Coeffs_nCoeffs))
    d1(:) = 0.0_R64
    d2(:) = 0.0_R64
    call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1, dCoeffs_=d1, h1_=h1)
    call Coeffs_ApplyH2FillRdm2(coeffs, rdm2_=rdm2, dCoeffs_=d2, h2_=h2)

  end subroutine

  !> Complex vector of length n with random real and imaginary parts.
  function RandomVector(n) result(res)
    integer(I32), intent(in) :: n
    complex(R64), allocatable :: res(:)

    real(R64) :: re(n), im(n)

    call random_number(re)
    call random_number(im)
    res = cmplx(re - 0.5_R64, im - 0.5_R64, kind=R64)

  end function

  !> Whether the total M (declared m times the occupation, summed over all
  !> orbitals) of each configuration is totalM.
  function InSector(totalM) result(res)
    integer(I32), intent(in) :: totalM
    logical, allocatable :: res(:)

    integer(I32) :: configurations(Method_Mb_nBodyTypes)
    integer(I32) :: iCoeff, bt, i, shift, m, iCNew
    real(R64) :: occupation

    allocate (res(Coeffs_nCoeffs))
    do iCoeff = 1, Coeffs_nCoeffs
      call Coeffs_ConfigurationsFromIndex(configurations, iCoeff)
      m = 0
      do bt = 1, Method_Mb_nBodyTypes
        shift = sum(Method_Mb_OrbBased_nOrbs(1:bt - 1))
        do i = 1, Method_Mb_OrbBased_nOrbs(bt)
          call configList(bt) % e % ExciteConfiguration(iCNew, occupation, [i], [i], configurations(bt))
          if (iCNew .ne. 0) m = m + nint(occupation) * Method_Mb_OrbBased_mOfOrb(shift + i)
        end do
      end do
      res(iCoeff) = m .eq. totalM
    end do

  end function

end program T_Coeffs_Generic_Sector
