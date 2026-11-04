!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Reader of the dipole file of a kick run (FemKick.inc, RunKick), shared by
! T_Fem_07_H2KickDouble and T_Fem_09_H2KickMeanFieldReference, which compare
! with the dipoleKick.dat of T_Fem_06_H2Kick (included after `contains`).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Header and dipole lines of a kick file of RunKick (FemKick.inc): the
!> header line and its run parameters, and the times and <z> (columns 1 and 4
!> of the lines of PrinterObservableFem_DumpDipole), indexed from 0. Stops if
!> the file or its header is missing (stale or foreign file).
subroutine ReadKickFile(header, kappaFile, dtFile, endTimeFile, energyGroundFile, nStepsFile, nPointsFile, &
                        fileTimes, fileDipoleZ, fileName)
  character(len=:), intent(out), allocatable :: header
  real(R64), intent(out) :: kappaFile
  real(R64), intent(out) :: dtFile
  real(R64), intent(out) :: endTimeFile
  real(R64), intent(out) :: energyGroundFile
  integer(I32), intent(out) :: nStepsFile
  integer(I32), intent(out) :: nPointsFile
  real(R64), intent(out), allocatable :: fileTimes(:)
  real(R64), intent(out), allocatable :: fileDipoleZ(:)
  character(len=*), intent(in) :: fileName

  character(len=1024) :: record
  real(R64) :: line(4)
  integer(I32) :: io, istat, nLines, i, equalsPos

  open (newunit=io, file=fileName, status="old", action="read", iostat=istat)
  if (istat .ne. 0) then
    write (*, '(1X, 3A)') "ReadKickFile: ", fileName, " missing"
    error stop "ReadKickFile: kick dipole file missing (run T_Fem_06_H2Kick first)"
  end if

  ! Header: "# T_Fem kick (<date>): ... = kappa timeStep endTime E_el(0) nSteps nPoints"
  read (io, '(A)', iostat=istat) record
  equalsPos = index(record, "=")
  if (istat .ne. 0 .or. index(record, "# T_Fem kick (") .ne. 1 .or. equalsPos .eq. 0) then
    write (*, '(1X, 3A)') "ReadKickFile: ", fileName, " has no kick header"
    error stop "ReadKickFile: kick dipole file without header (stale file; run T_Fem_06_H2Kick first)"
  end if
  read (record(equalsPos + 1:), *, iostat=istat) kappaFile, dtFile, endTimeFile, energyGroundFile, nStepsFile, nPointsFile
  if (istat .ne. 0) error stop "ReadKickFile: unreadable kick header"
  header = trim(record)

  nLines = 0
  do
    read (io, *, iostat=istat) line
    if (istat .ne. 0) exit
    nLines = nLines + 1
  end do
  rewind (io)
  read (io, '(A)') record
  allocate (fileTimes(0:nLines - 1), fileDipoleZ(0:nLines - 1))
  do i = 0, nLines - 1
    read (io, *) line
    fileTimes(i) = line(1)
    fileDipoleZ(i) = line(4)
  end do
  close (io)

end subroutine
