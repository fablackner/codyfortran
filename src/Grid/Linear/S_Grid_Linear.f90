! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Submodule implementing a linear grid in one dimension.
submodule(M_Grid_Linear) S_Grid_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Sets up the linear grid based on parameters from the JSON configuration
  module subroutine Grid_Linear_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Linear_Const
    use M_Grid_Linear_Fedvr

    call Say_Fabricate("grid.linear")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Get configuration parameters from JSON
    Grid_Linear_xmin = Json_Get("grid.linear.xmin", -10.0_R64)
    Grid_Linear_xmax = Json_Get("grid.linear.xmax", 10.0_R64)

    ! Set procedure pointers
    Grid_InnerProduct => InnerProduct

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("grid.linear.const")) then
      call Grid_Linear_Const_Fabricate

    else if (Json_GetExistence("grid.linear.fedvr")) then
      call Grid_Linear_Fedvr_Fabricate

    else
      error stop "grid.linear is missing one of: const, fedvr"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Computes the inner product of two fields
  pure function InnerProduct(fConjg, f) result(res)

    complex(R64), intent(in), contiguous, dimension(:) :: fConjg
    complex(R64), intent(in), contiguous, dimension(:) :: f
    complex(R64) :: res

    res = sum(conjg(fConjg(:)) * f(:) * Grid_Linear_weights(:))

  end function

end submodule
