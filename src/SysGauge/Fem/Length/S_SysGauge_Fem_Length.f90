! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysGauge_Fem_Length) S_SysGauge_Fem_Length

  implicit none

  !-----------------------------------------------------------------------------
  ! local data
  !-----------------------------------------------------------------------------

  !> e_pol . r at the nodes (Grid_nPoints), built in Setup
  !> (Grid_Fem_FillNodeCoordinate).
  real(R64), allocatable :: polarizationDotR(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysGauge_Fem_Length_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysGauge

    real(R64), allocatable :: polarization(:)
    real(R64) :: norm

    call Say_Fabricate("sysGauge.fem.length")

    !------------------------------------
    ! pulse parameters
    !------------------------------------

    SysGauge_Fem_Length_fieldStrength = Json_Get("sysGauge.fem.length.fieldStrength", 0.01_R64)
    SysGauge_Fem_Length_omega = Json_Get("sysGauge.fem.length.omega", 0.057_R64)
    SysGauge_Fem_Length_nCycles = Json_Get("sysGauge.fem.length.nCycles", 3.0_R64)
    SysGauge_Fem_Length_cep = Json_Get("sysGauge.fem.length.cep", 0.0_R64)
    SysGauge_Fem_Length_tStart = Json_Get("sysGauge.fem.length.tStart", 0.0_R64)
    polarization = Json_Get("sysGauge.fem.length.polarization", [0.0_R64, 0.0_R64, 1.0_R64])

    if (.not. (SysGauge_Fem_Length_omega > 0.0_R64)) error stop "sysGauge.fem.length.omega must be > 0"
    if (.not. (SysGauge_Fem_Length_nCycles > 0.0_R64)) error stop "sysGauge.fem.length.nCycles must be > 0"

    if (size(polarization) .ne. 3) error stop "sysGauge.fem.length.polarization needs three components"
    norm = norm2(polarization)
    if (.not. (norm > 0.0_R64)) error stop "sysGauge.fem.length.polarization must be nonzero"
    SysGauge_Fem_Length_polarization(:) = polarization(:) / norm

    !------------------------------------
    ! flags and procedure pointers
    !------------------------------------

    SysGauge_timeIndependentQ = .false.
    SysGauge_bodyTypeIndependentQ = .true.

    SysGauge_Setup => Setup
    SysGauge_MultiplyWithGaugeOp => MultiplyWithGaugeOp

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Electric field E(t) of the sin^2-envelope pulse (zero outside the pulse),
  !> the pulse of sysPotential.ylm.coulombLaser (Pulse_Sin2Field).
  module function SysGauge_Fem_Length_FieldAmplitude(time) result(res)
    use M_Utils_Pulse

    real(R64), intent(in) :: time
    real(R64) :: res

    res = Pulse_Sin2Field(time, SysGauge_Fem_Length_fieldStrength, SysGauge_Fem_Length_omega, &
                          SysGauge_Fem_Length_nCycles, SysGauge_Fem_Length_cep, SysGauge_Fem_Length_tStart)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Cache e_pol . r at the nodes (after Grid_Setup).
  subroutine Setup
    use M_Utils_Say
    use M_Grid_Fem

    call Say_Setup("sysGauge.fem.length")

    call Grid_Fem_FillNodeCoordinate(polarizationDotR, SysGauge_Fem_Length_polarization)

    write (*, '(A, 3F10.6)') "  sysGauge.fem.length: polarization = ", SysGauge_Fem_Length_polarization(:)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> dOrb = E(t) (e_pol . r) orb, point-wise (zero outside the pulse); the
  !> coupling is the same for every body type (bt_ is ignored). Reentrant.
  subroutine MultiplyWithGaugeOp(dOrb, orb, time, bt_)
    !$ use omp_lib, only: omp_in_parallel
    use M_Utils_UnusedVariables

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous  :: orb(:)
    real(R64), intent(in)                 :: time
    integer(I32), intent(in), optional    :: bt_

    real(R64) :: field
    integer(I32) :: i

    if (.false.) call UnusedVariables_Mark(bt_)

    if (.not. allocated(polarizationDotR)) then
      error stop "sysGauge.fem.length: call SysGauge_Setup before SysGauge_MultiplyWithGaugeOp"
    end if

    field = SysGauge_Fem_Length_FieldAmplitude(time)

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
