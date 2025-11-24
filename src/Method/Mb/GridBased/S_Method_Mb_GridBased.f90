! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_GridBased) S_Method_Mb_GridBased

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_GridBased_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Method
    use M_Method_Mb
    use M_Method_Mb_GridBased_Full

    call Say_Fabricate("method.gridBased")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Method_GetEnergy => GetEnergy
    Method_Mb_GridBased_ApplyKineticOp => ApplyKineticOp
    Method_Mb_GridBased_ApplyPotentialOp => ApplyPotentialOp
    Method_Mb_GridBased_ApplyInteractionOp => ApplyInteractionOp

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("method.mb.gridBased.full")) then
      call Method_Mb_GridBased_Full_Fabricate

    else
      error stop "method.mb.gridBased is missing one of: full"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function GetEnergy(time) result(res)
    use M_Grid
    use M_Method

    real(R64)                :: res
    real(R64), intent(in)    :: time

    complex(R64), allocatable :: dState(:)

    allocate (dState, mold=Method_state)

    call Method_TimeDerivative(dState, Method_state, time)

    res = real(Grid_InnerProduct(dState, Method_state), kind=R64)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyKineticOp(dSlice, slice, time, bt_)
    use M_Grid
    use M_Method
    use M_Method_Mb_GridBased
    use M_SysKinetic

    complex(R64), intent(inout), contiguous :: dSlice(:)
    complex(R64), intent(in), contiguous    :: slice(:)
    real(R64), intent(in)                   :: time
    integer(I32), intent(in), optional      :: bt_

    complex(R64), allocatable :: sliceTmp(:)

    allocate (sliceTmp(size(slice)))

    sliceTmp(:) = 0.0_R64
    call SysKinetic_MultiplyWithKineticOp(sliceTmp, slice(:), time, bt_)

    dSlice(:) = dSlice(:) + sliceTmp(:)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyPotentialOp(dSlice, slice, time, bt_)
    use M_Grid
    use M_SysPotential

    complex(R64), intent(inout), contiguous :: dSlice(:)
    complex(R64), intent(in), contiguous    :: slice(:)
    real(R64), intent(in)                   :: time
    integer(I32), intent(in), optional      :: bt_

    complex(R64), allocatable :: externalPotential(:)
    complex(R64), allocatable :: orbTmp(:)
    integer(I32) :: nG

    nG = Grid_nPoints
    allocate (orbTmp(nG))

    call SysPotential_FillExternalPotential(externalPotential, time, bt_)
    call SysPotential_MultiplyWithExternalPotential(orbTmp, externalPotential, slice)

    dSlice(:) = dSlice(:) + orbTmp(:)

    deallocate (orbTmp)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyInteractionOp(dSlice2d, slice2d, bt1, bt2, time)
    use M_Grid
    use M_SysInteraction

    complex(R64), intent(inout), contiguous :: dSlice2d(:, :)
    complex(R64), intent(in), contiguous :: slice2d(:, :)
    integer(I32), intent(in) :: bt1
    integer(I32), intent(in) :: bt2
    real(R64), intent(in) :: time

    integer(I32) :: nG, iGrid2
    complex(R64), allocatable :: src(:), interactionPotential(:), orbTmp(:)

    nG = Grid_nPoints

    allocate (orbTmp(nG))
    allocate (src(nG))

    do iGrid2 = 1, nG
      src = (0.0_R64, 0.0_R64)
      src(iGrid2) = (1.0_R64, 0.0_R64)
      call SysInteraction_FillInteractionPotential(interactionPotential, src, time, bt1, bt2)
      call SysInteraction_MultiplyWithInteractionPotential(orbTmp, interactionPotential, slice2d(:, iGrid2))
      dSlice2d(:, iGrid2) = dSlice2d(:, iGrid2) + orbTmp(:)
    end do

    deallocate (src, orbTmp, interactionPotential)

  end subroutine

end submodule
