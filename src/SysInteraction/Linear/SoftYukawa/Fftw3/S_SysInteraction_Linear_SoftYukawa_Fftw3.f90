! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief FFT-based convolution implementation of SoftYukawa.
submodule(M_SysInteraction_Linear_SoftYukawa_Fftw) S_SysInteraction_Linear_SoftYukawa_Fftw
  use M_Utils_ConvolutionFftw

  implicit none

  !> Convolution context holding FFT plans and transformed kernel
  type(T_ConvolutionFftw_Ctx) :: convolutionCtx

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind the FFT-based SoftYukawa implementation.
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
  !> @brief Initialize FFT plans and transform the kernel.
  !>
  !> Pre-transforms the SoftYukawa kernel to Fourier space for fast convolution.
  subroutine Setup
    use M_Utils_Constants
    use M_SysInteraction_Linear_SoftYukawa
    use M_Grid
    use M_Grid_Linear
    implicit none

    integer(I32) :: nG
    real(R64)    :: dx

    nG = Grid_nPoints
    dx = Grid_Linear_xCoord(2) - Grid_Linear_xCoord(1)

    call ConvolutionFftw_CreateCtx(convolutionCtx, nG, dx, SysInteraction_Linear_SoftYukawa_Interaction)

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Compute interaction potential via FFT convolution.
  !>
  !> Evaluates V = IFFT( K̃ * FFT(ρ) ) where K̃ is the pre-transformed kernel.
  !>
  !> @param[out] interactionPotential  Resulting potential V(x)
  !> @param[in]  src                   Weighted source density ρ(x)·w(x)
  !> @param[in]  time                  Physical time (unused)
  !> @param[in]  bt1_                  Target body type (unused)
  !> @param[in]  bt2_                  Source body type (unused)
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

    call ConvolutionFftw_Apply(interactionPotential, src, convolutionCtx)

  end subroutine

end submodule
