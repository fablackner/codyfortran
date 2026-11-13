! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysGauge_Length) S_SysGauge_Length

  implicit none

  !-----------------------------------------------------------------------------
  ! local data
  !-----------------------------------------------------------------------------

  !> JSON path of the model (messages).
  character(len=:), allocatable :: keyPath
  !> Coordinate filler of the grid family.
  procedure(I_SysGauge_Length_FillCoordinate), pointer :: FillPolarizationDotR => null()
  !> e_pol . r at the nodes (Grid_nPoints), built in Setup.
  real(R64), allocatable :: polarizationDotR(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysGauge_Length_Fabricate(path, FillCoordinate)
    use M_Utils_Json
    use M_Utils_Say
    use M_SysGauge

    character(len=*), intent(in) :: path
    procedure(I_SysGauge_Length_FillCoordinate) :: FillCoordinate

    real(R64), allocatable :: polarization(:)

    call Say_Fabricate(path)

    keyPath = path
    SysGauge_Length_fieldStrength = Json_Get(path//".fieldStrength", 0.01_R64)
    SysGauge_Length_omega = Json_Get(path//".omega", 0.057_R64)
    SysGauge_Length_nCycles = Json_Get(path//".nCycles", 3.0_R64)
    SysGauge_Length_cep = Json_Get(path//".cep", 0.0_R64)
    SysGauge_Length_tStart = Json_Get(path//".tStart", 0.0_R64)
    polarization = Json_Get(path//".polarization", [0.0_R64, 0.0_R64, 1.0_R64])

    if (.not. (SysGauge_Length_omega > 0.0_R64)) error stop path//".omega must be > 0"
    if (.not. (SysGauge_Length_nCycles > 0.0_R64)) error stop path//".nCycles must be > 0"
    if (size(polarization) .ne. 3) error stop path//".polarization needs three components"
    if (.not. (norm2(polarization) > 0.0_R64)) error stop path//".polarization must be nonzero"
    SysGauge_Length_polarization(:) = polarization(:) / norm2(polarization)

    if (allocated(polarizationDotR)) deallocate (polarizationDotR)
    FillPolarizationDotR => FillCoordinate

    SysGauge_timeIndependentQ = .false.
    SysGauge_bodyTypeIndependentQ = .true.

    SysGauge_Setup => Setup
    SysGauge_MultiplyWithGaugeOp => MultiplyWithGaugeOp

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function SysGauge_Length_FieldAmplitude(time) result(res)
    use M_Utils_Pulse

    real(R64), intent(in) :: time
    real(R64) :: res

    res = Pulse_Sin2Field(time, SysGauge_Length_fieldStrength, SysGauge_Length_omega, &
                          SysGauge_Length_nCycles, SysGauge_Length_cep, SysGauge_Length_tStart)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Cache e_pol . r at the nodes (after Grid_Setup).
  subroutine Setup
    use M_Utils_Say

    call Say_Setup(keyPath)

    call FillPolarizationDotR(polarizationDotR, SysGauge_Length_polarization)

    write (*, '(A, 3F10.6)') "  "//keyPath//": polarization = ", SysGauge_Length_polarization(:)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> dOrb = E(t) (e_pol . r) orb, point-wise (zero outside the pulse); the
  !> coupling is the same for every body type (bt_ is ignored). Reentrant.
  subroutine MultiplyWithGaugeOp(dOrb, orb, time, bt_)
!$  use omp_lib, only: omp_in_parallel
    use M_Utils_UnusedVariables

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous  :: orb(:)
    real(R64), intent(in)                 :: time
    integer(I32), intent(in), optional    :: bt_

    real(R64) :: field
    integer(I32) :: i

    if (.false.) call UnusedVariables_Mark(bt_)

    if (.not. allocated(polarizationDotR)) error stop keyPath//": call SysGauge_Setup before SysGauge_MultiplyWithGaugeOp"

    field = SysGauge_Length_FieldAmplitude(time)

    if (abs(field) <= 0.0_R64) then
      dOrb(:) = (0.0_R64, 0.0_R64)
      return
    end if

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(dOrb)
      dOrb(i) = (field * polarizationDotR(i)) * orb(i)
    end do
    !$omp end parallel do

  end subroutine

end submodule
