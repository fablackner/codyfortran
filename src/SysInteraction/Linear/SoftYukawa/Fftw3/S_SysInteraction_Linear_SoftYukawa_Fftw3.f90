! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysInteraction_Linear_SoftYukawa_Fftw) S_SysInteraction_Linear_SoftYukawa_Fftw
  use M_Utils_ConvolutionFftw

  implicit none

  type(T_ConvolutionFftw_Ctx) :: convolutionCtx

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_Linear_SoftYukawa_Fftw_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Linear
    use M_SysInteraction_Linear_SoftYukawa

    call Say_Fabricate("sysInteraction.linear.softYukawa.fftw")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_FillInteractionPotential => FillInteractionPotential
    SysInteraction_Setup => Setup

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Constants
    use M_SysInteraction_Linear_SoftYukawa
    use M_Grid
    use M_Grid_Linear
    implicit none

    integer(I32) :: nG
    real(R64)    :: dx

    ! --- Obtain grid info ---
    nG = Grid_nPoints
    dx = Grid_Linear_xCoord(2) - Grid_Linear_xCoord(1)

    ! --- Initialize convolution with the SoftYukawa kernel function ---
    call ConvolutionFftw_CreateCtx(convolutionCtx, nG, dx, SysInteraction_Linear_SoftYukawa_Interaction)

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillInteractionPotential(interactionPotential, src, time, bt1_, bt2_)
    use M_Utils_UnusedVariables
    use M_Utils_Constants
    use M_Grid
    use M_Utils_ConvolutionFftw
    implicit none

    complex(R64), intent(out), allocatable :: interactionPotential(:)
    complex(R64), intent(in), contiguous :: src(:)
    real(R64), intent(in)  :: time
    integer(I32), intent(in), optional  :: bt1_
    integer(I32), intent(in), optional  :: bt2_

    if (.false.) call UnusedVariables_Mark(time, bt1_, bt2_)

    if (.not. allocated(interactionPotential)) allocate (interactionPotential(Grid_nPoints))

    ! --- Apply convolution utility ---
    call ConvolutionFftw_Apply(interactionPotential, src, convolutionCtx)

  end subroutine

end submodule
