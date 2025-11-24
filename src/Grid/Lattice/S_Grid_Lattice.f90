! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Lattice) S_Grid_Lattice

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Lattice_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid

    call Say_Fabricate("grid.lattice")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Grid_Lattice_xSize = Json_Get("grid.lattice.xSize", 1)
    Grid_Lattice_ySize = Json_Get("grid.lattice.ySize", 1)
    Grid_Lattice_zSize = Json_Get("grid.lattice.zSize", 1)

    Grid_nPoints = Grid_Lattice_xSize * Grid_Lattice_ySize * Grid_Lattice_zSize

    Grid_Lattice_xPeriodicQ = Json_Get("grid.lattice.xPeriodicQ", .false.)
    Grid_Lattice_yPeriodicQ = Json_Get("grid.lattice.yPeriodicQ", .false.)
    Grid_Lattice_zPeriodicQ = Json_Get("grid.lattice.zPeriodicQ", .false.)

    Grid_Setup => Setup
    Grid_InnerProduct => InnerProduct

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid

    integer(I32) :: iGrid
    integer(I32) :: ix, iy, iz

    call Say_Setup("grid.lattice")

    allocate (Grid_Lattice_code(Grid_Lattice_xSize, Grid_Lattice_ySize, Grid_Lattice_zSize))

    iGrid = 0
    do iz = 1, Grid_Lattice_zSize
      do iy = 1, Grid_Lattice_ySize
        do ix = 1, Grid_Lattice_xSize
          iGrid = iGrid + 1

          Grid_Lattice_code(ix, iy, iz) = iGrid

        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Computes the inner product of two fields
  pure function InnerProduct(fConjg, f) result(res)

    complex(R64), intent(in), contiguous, dimension(:) :: fConjg
    complex(R64), intent(in), contiguous, dimension(:) :: f
    complex(R64) :: res

    !res = sum(conjg(fConjg(:)) * f(:))
    ! does the same as above, but (maybe) faster
    res = dot_product(fConjg(:), f(:))

  end function

end submodule
