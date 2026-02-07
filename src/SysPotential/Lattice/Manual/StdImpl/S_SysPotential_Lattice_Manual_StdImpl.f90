! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Lattice_Manual_StdImpl) S_SysPotential_Lattice_Manual_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Lattice_Manual_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential, only: SysPotential_FillExternalPotential, SysPotential_Setup

    call Say_Fabricate("sysPotential.lattice.manual.stdImpl")

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
    use M_SysPotential_Lattice_Manual

    integer(I32) :: nG, i, iSite

    call Say_Setup("sysPotential.lattice.manual.stdImpl")

    nG = Grid_nPoints

    allocate (SysPotential_Lattice_Manual_manualValues(nG))
    SysPotential_Lattice_Manual_manualValues = 0.0_R64

    if (allocated(SysPotential_Lattice_Manual_sites)) then
      if (size(SysPotential_Lattice_Manual_values) /= size(SysPotential_Lattice_Manual_sites)) then
        error stop "manual: values and sites arrays must have the same size"
      end if

      do i = 1, size(SysPotential_Lattice_Manual_sites)
        iSite = SysPotential_Lattice_Manual_sites(i)
        if (iSite < 1 .or. iSite > nG) error stop "manual: site index out of bounds"
        SysPotential_Lattice_Manual_manualValues(iSite) = SysPotential_Lattice_Manual_values(i)
      end do
    else
      if (size(SysPotential_Lattice_Manual_values) /= nG) then
        error stop "manual: values array size must match grid size when sites are not specified"
      end if

      SysPotential_Lattice_Manual_manualValues = SysPotential_Lattice_Manual_values
    end if

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillExternalPotential(externalPotential, time, bt_)
    use M_Utils_UnusedVariables
    use M_SysPotential_Lattice_Manual
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
          externalPotential(iGrid) = SysPotential_Lattice_Manual_manualValues(iGrid)
        end do
      end do
    end do

  end subroutine

end submodule
