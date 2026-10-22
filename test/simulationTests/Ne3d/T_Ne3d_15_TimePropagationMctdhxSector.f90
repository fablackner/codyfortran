!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> T_Ne3d_13_TimePropagationMctdhxDeclaredM for Ne6+ (Be-like, 2 + 2
!> electrons) in the same 5 restricted orbitals with all configurations per
!> spin (100 in all), so that the CI sector of total M = 0
!> (coeffs.generic.totalM, 34 of 100 configurations) is a proper subspace:
!> the MCSCF ground solver projects its CI eigenvectors on it, and the
!> Hamiltonian acts on the sector only. The 2s-2p near-degeneracy gives
!> natural occupations of 2e-2 and more, so the propagation does not amplify
!> round-off (a weakly occupied orbital, e.g. 3s for Ne, would).
!>
!> Checks:
!> - the coefficients outside the sector are exactly zero after every output
!>   step (h1, h2 and the gauge coupling conserve M)
!> Writes the ground-state energy and the dipole history dipoleSector.dat,
!> which T_Ne3d_16_TimePropagationMctdhxFullCi compares with the full space.
program T_Ne3d_15_TimePropagationMctdhxSector
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Coeffs
  use M_Method
  use testdrive, only: check, error_type

  implicit none

  real(R64)    :: time, energy, dipole
  integer(I32) :: nOtherM, io
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupMctdhxLaser("test/simulationTests/Ne3d/T_Ne3d_15_TimePropagationMctdhxSector.json")

  !==========================================
  call Say_Section("converge mcscf ground state")
  !==========================================

  energy = ConvergedGroundStateEnergy()

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  nOtherM = OtherMCoeffs(Method_state(1:Coeffs_nCoeffs))

  open (newunit=io, file="dipoleSector.dat", status="replace")
  write (io, '(ES26.17)') energy

  time = Json_Get("program.startTime", 0.0_R64)
  do while (BeforeEndTimeQ(time))
    call PropagateOutputStep(time, dipole, "dipoleVelocity.dat")
    nOtherM = nOtherM + OtherMCoeffs(Method_state(1:Coeffs_nCoeffs))
    write (io, '(2ES26.17)') time, dipole
  end do

  close (io)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "final dipole: ", dipole
  print *, "nonzero coefficients outside the sector: ", nOtherM

  ! Exact zeros, not round-off
  call check(error, nOtherM, 0)
  if (allocated(error)) error stop "T_Ne3d_15_TimePropagationMctdhxSector failure: coefficients left the sector"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'Ne3dMctdhxLaser.inc'

  !> Number of nonzero CI coefficients of configurations of total M other
  !> than 0 (declared m over the occupied orbitals, fermionic bit codes).
  function OtherMCoeffs(coeffs) result(res)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_ConfigList

    complex(R64), intent(in), contiguous :: coeffs(:)
    integer(I32) :: res

    integer(I32) :: configurations(Method_Mb_nBodyTypes)
    integer(I32) :: iCoeff, bt, i, shift, totalM

    res = 0
    do iCoeff = 1, size(coeffs)
      call Coeffs_ConfigurationsFromIndex(configurations, iCoeff)
      totalM = 0
      do bt = 1, Method_Mb_nBodyTypes
        shift = sum(Method_Mb_OrbBased_nOrbs(1:bt - 1))
        do i = 1, Method_Mb_OrbBased_nOrbs(bt)
          if (btest(configList(bt) % e % codeFromConfig(configurations(bt)), i - 1)) &
            totalM = totalM + Method_Mb_OrbBased_mOfOrb(shift + i)
        end do
      end do
      if (totalM .ne. 0 .and. coeffs(iCoeff) .ne. (0.0_R64, 0.0_R64)) res = res + 1
    end do

  end function

end program
