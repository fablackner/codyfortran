! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Lattice_RandomGauss_StdImpl) S_SysPotential_Lattice_RandomGauss_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Lattice_RandomGauss_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential, only: SysPotential_FillExternalPotential, SysPotential_Setup

    call Say_Fabricate("sysPotential.lattice.randomGauss.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_FillExternalPotential => FillExternalPotential
    SysPotential_Setup => Setup

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_SysPotential_Lattice_RandomGauss

    real(R64) :: meanV, stdV
    integer(I32) :: seedValue
    real(R64) :: u1, u2, r, theta, pi
    integer(I32) :: nG, i
    integer :: seedSize, seedClock
    integer, allocatable :: seed(:)

    call Say_Setup("sysPotential.lattice.randomGauss.stdImpl")

    ! Pull config
    meanV = SysPotential_Lattice_RandomGauss_meanValue
    stdV = SysPotential_Lattice_RandomGauss_stdValue
    seedValue = SysPotential_Lattice_RandomGauss_seed
    nG = Grid_nPoints

    ! RNG seeding
    call random_seed(size=seedSize)
    allocate (seed(seedSize))
    if (seedValue .eq. -1) then
      call system_clock(count=seedClock)
      seed = seedClock + [(i, i=1, seedSize)] * 37
    else
      seed = seedValue
    end if
    call random_seed(put=seed)
    deallocate (seed)

    allocate (SysPotential_Lattice_RandomGauss_randomGaussValues(nG))

    if (stdV <= 0.0_R64) error stop 'randomGauss: stdValue must be > 0 for gaussian distribution'
    pi = acos(-1.0_R64)
    i = 1
    do while (i <= nG)
      call random_number(u1)
      call random_number(u2)
      if (u1 <= 1.0e-12_R64) cycle
      r = sqrt(-2.0_R64 * log(u1))
      theta = 2.0_R64 * pi * u2
      SysPotential_Lattice_RandomGauss_randomGaussValues(i) = meanV + stdV * (r * cos(theta))
      if (i + 1 <= nG) then
        SysPotential_Lattice_RandomGauss_randomGaussValues(i + 1) = meanV + stdV * (r * sin(theta))
      end if
      i = i + 2
    end do

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillExternalPotential(externalPotential, time, bt_)
    use M_Utils_UnusedVariables
    use M_SysPotential_Lattice_RandomGauss
    use M_Grid
    use M_Grid_Lattice

    complex(R64), intent(out), allocatable :: externalPotential(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: ix, iy, iz, iGrid

    if (.false.) call UnusedVariables_Mark(time, bt_)

    if (.not. allocated(externalPotential)) allocate (externalPotential(Grid_nPoints))

    ! write cached values into externalPotential in lattice order
    do iz = 1, Grid_Lattice_zSize
      do iy = 1, Grid_Lattice_ySize
        do ix = 1, Grid_Lattice_xSize
          iGrid = Grid_Lattice_code(ix, iy, iz)
          externalPotential(iGrid) = SysPotential_Lattice_RandomGauss_randomGaussValues(iGrid)
        end do
      end do
    end do

  end subroutine

end submodule
