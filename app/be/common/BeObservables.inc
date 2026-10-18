!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Observables of the beryllium laser programs (included after `contains`).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Append one line to observables.dat and natocc.dat (replaced if newFilesQ):
!>   observables.dat: time, E(t), A(t), <z>, -Z<z/r^3>, electrons in the box,
!>                    energy
!>   natocc.dat:      time, natural occupations of the spin-up 1-RDM
!> Absorbed orbitals are no longer orthonormal: the electron number is
!> Tr(rho S) and the occupations are the eigenvalues of S^1/2 rho S^1/2 with
!> the orbital overlap matrix S (i.e. of the density matrix inside the box).
subroutine WriteLaserObservables(rdm1, time, charge, energy, newFilesQ)
  use M_Utils_LapackLib
  use M_Utils_PrinterObservableYlm
  use M_Grid
  use M_Orbs
  use M_SysGauge_Ylm_Velocity

  complex(R64), intent(in), contiguous :: rdm1(:, :)
  real(R64), intent(in) :: time
  real(R64), intent(in) :: charge
  real(R64), intent(in) :: energy
  logical, intent(in) :: newFilesQ

  real(R64) :: dipole, acceleration, nElectrons
  real(R64), allocatable :: evals(:), occ(:)
  complex(R64), allocatable :: s(:, :), evecs(:, :), sSqrt(:, :), rhoUp(:, :)
  integer(I32) :: nOS, i, j, io
  character(len=7) :: openStatus

  nOS = Orbs_nOrbsInState

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
  sSqrt = matmul(evecs, spread(sqrt(max(evals, 0.0_R64)), 2, nOS) * conjg(transpose(evecs)))
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
