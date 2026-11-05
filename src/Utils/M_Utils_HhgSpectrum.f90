! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> High-harmonic spectrum of a dipole time series (stateless).
!>
!> With the Hann window w(t) = sin^2(pi (t - t0)/T) over the whole series
!> t0 .. t0 + T and harmonic orders q of the carrier omega0, omega = q omega0:
!>   acceleration form: |FT[a(t) w(t)]|^2
!>   length form:       omega^4 |FT[(z(t) - <z>) w(t)]|^2
!> with FT[f](omega) = sum_k f(t_k) exp(-i omega t_k) dt on an equidistant
!> grid. Orders run from dOrder to maxOrder in steps of dOrder and stop at
!> the Nyquist frequency pi/dt.
module M_Utils_HhgSpectrum
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Read the columns colTime, colDipole and (optional) colAcceleration_ of
  !> the whitespace-separated table inFile (lines starting with # skipped) and
  !> write the spectrum to outFile (HhgSpectrum_Write).
  subroutine HhgSpectrum_WriteFromTable(outFile, inFile, omega0, colTime, colDipole, &
                                        colAcceleration_, maxOrder_, dOrder_)
    character(len=*), intent(in) :: outFile
    character(len=*), intent(in) :: inFile
    real(R64), intent(in) :: omega0
    integer(I32), intent(in) :: colTime
    integer(I32), intent(in) :: colDipole
    integer(I32), intent(in), optional :: colAcceleration_
    real(R64), intent(in), optional :: maxOrder_
    real(R64), intent(in), optional :: dOrder_

    real(R64), allocatable :: time(:), dipole(:), acceleration(:), row(:)
    character(len=4096) :: line
    integer(I32) :: io, ios, nCols, nRows

    nCols = max(colTime, colDipole)
    if (present(colAcceleration_)) nCols = max(nCols, colAcceleration_)
    allocate (row(nCols), time(0), dipole(0), acceleration(0))

    open (newunit=io, file=inFile, status="old", action="read")
    do
      read (io, '(A)', iostat=ios) line
      if (ios .ne. 0) exit
      if (len_trim(line) .eq. 0) cycle
      if (line(1:1) .eq. "#") cycle
      read (line, *) row
      time = [time, row(colTime)]
      dipole = [dipole, row(colDipole)]
      if (present(colAcceleration_)) acceleration = [acceleration, row(colAcceleration_)]
    end do
    close (io)

    nRows = size(time)
    if (nRows < 2) error stop "HhgSpectrum_WriteFromTable: fewer than two rows"

    if (present(colAcceleration_)) then
      call HhgSpectrum_Write(outFile, time, dipole, omega0, acceleration_=acceleration, &
                             maxOrder_=maxOrder_, dOrder_=dOrder_)
    else
      call HhgSpectrum_Write(outFile, time, dipole, omega0, maxOrder_=maxOrder_, dOrder_=dOrder_)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Write the spectrum of the equidistant series time, dipole <z> and
  !> (optional) dipole acceleration to outFile: harmonic order, then the
  !> acceleration form (if given) and the length form. maxOrder_ defaults
  !> to 300, dOrder_ to 0.1.
  subroutine HhgSpectrum_Write(outFile, time, dipole, omega0, acceleration_, maxOrder_, dOrder_)
    use M_Utils_Constants, only: PI

    character(len=*), intent(in) :: outFile
    real(R64), intent(in), contiguous :: time(:)
    real(R64), intent(in), contiguous :: dipole(:)
    real(R64), intent(in) :: omega0
    real(R64), intent(in), contiguous, optional :: acceleration_(:)
    real(R64), intent(in), optional :: maxOrder_
    real(R64), intent(in), optional :: dOrder_

    real(R64) :: maxOrder, dOrder, span, dt, dipoleMean, omega
    real(R64), allocatable :: window(:)
    complex(R64), allocatable :: phase(:)
    complex(R64) :: ftDipole, ftAcceleration
    integer(I32) :: n, iOrder, io

    maxOrder = 300.0_R64
    if (present(maxOrder_)) maxOrder = maxOrder_
    dOrder = 0.1_R64
    if (present(dOrder_)) dOrder = dOrder_

    n = size(time)
    span = time(n) - time(1)
    dt = span / (n - 1)
    window = sin(PI * (time - time(1)) / span)**2
    dipoleMean = sum(dipole) / n

    open (newunit=io, file=outFile, status="replace", action="write")
    if (present(acceleration_)) then
      write (io, '(A)') "# harmonic order, |a(w)|^2, w^4 |z(w)|^2"
    else
      write (io, '(A)') "# harmonic order, w^4 |z(w)|^2"
    end if

    iOrder = 1
    do while (iOrder * dOrder <= maxOrder * (1.0_R64 + 1e-12_R64))
      omega = iOrder * dOrder * omega0
      if (omega >= PI / dt) exit

      phase = exp(cmplx(0.0_R64, -omega * time, kind=R64)) * window * dt
      ftDipole = sum((dipole - dipoleMean) * phase)

      if (present(acceleration_)) then
        ftAcceleration = sum(acceleration_ * phase)
        write (io, '(F10.2, 2ES16.6E3)') iOrder * dOrder, abs(ftAcceleration)**2, omega**4 * abs(ftDipole)**2
      else
        write (io, '(F10.2, ES16.6E3)') iOrder * dOrder, omega**4 * abs(ftDipole)**2
      end if

      iOrder = iOrder + 1
    end do

    close (io)

  end subroutine

end module
