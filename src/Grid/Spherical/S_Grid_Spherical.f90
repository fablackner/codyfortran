! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Spherical) S_Grid_Spherical

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Spherical_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Spherical_Const

    call Say_Fabricate("grid.spherical")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Read radial grid parameters from JSON
    Grid_Spherical_rmin = Json_Get("grid.spatial.constRadial3d.rmin", 0.0_R64)
    Grid_Spherical_rmax = Json_Get("grid.spatial.constRadial3d.rmax", 20.0_R64)

    Grid_InnerProduct => InnerProduct

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("grid.spherical.const")) then
      call Grid_Spherical_Const_Fabricate
    else
      error stop "grid.spherical is missing one of: const"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Computes the inner product of two fields
  pure function InnerProduct(fConjg, f) result(res)

    complex(R64), intent(in), contiguous, dimension(:) :: fConjg
    complex(R64), intent(in), contiguous, dimension(:) :: f
    complex(R64) :: res

    res = sum(conjg(fConjg(:)) * f(:) * Grid_Spherical_weights(:))

  end function

end submodule
