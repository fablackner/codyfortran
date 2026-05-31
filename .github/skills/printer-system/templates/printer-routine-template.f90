! Template for a new printer routine in M_Utils_Printer*.f90 modules
subroutine PrinterX_DumpObservable(data, filename, toScreenQ, time)
  use M_Utils_Types

  complex(R64), intent(in), contiguous :: data(:, :)
  character(len=*), intent(in) :: filename
  logical, intent(in) :: toScreenQ
  real(R64), intent(in) :: time

  integer(I32) :: io, istat
  character(256) :: msg
  complex(R64) :: value
  logical :: needOutput

  needOutput = toScreenQ .or. (filename .ne. "")
  if (.not. needOutput) return

  ! Compute observable from data here.
  value = 0.0_R64

  ! Print to Screen
  if (toScreenQ) then
    write (*, '(A, 2E20.10E3)') "Observable: ", value
  end if

  ! Print to file
  if (filename .eq. "") return

  open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
  if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
  if (istat .ne. 0) error stop

  write (io, '(F10.6, 2E20.10E3)') time, value

  close (io)

end subroutine
