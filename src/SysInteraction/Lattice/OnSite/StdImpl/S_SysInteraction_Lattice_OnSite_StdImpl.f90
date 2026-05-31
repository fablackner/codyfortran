! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Standard implementation of the on-site lattice interaction.
submodule(M_SysInteraction_Lattice_OnSite_StdImpl) S_SysInteraction_Lattice_OnSite_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind the standard on-site implementation.
  module subroutine SysInteraction_Lattice_OnSite_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction

    call Say_Fabricate("sysInteraction.lattice.onSite.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_FillInteractionPotential => FillInteractionPotential

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Compute interaction potential from source density.
  !>
  !> For on-site interaction: V(i) = U * ρ(i)
  !> This is a simple scaling—no convolution or Poisson solve needed.
  !>
  !> @param[out] interactionPotential  Resulting potential V(i)
  !> @param[in]  src                   Source density ρ(i)
  !> @param[in]  time                  Physical time (unused, time-independent)
  !> @param[in]  bt1_                  Target body type (unused)
  !> @param[in]  bt2_                  Source body type (unused)
  subroutine FillInteractionPotential(interactionPotential, src, time, bt1_, bt2_)
    use M_Utils_UnusedVariables
    use M_SysInteraction_Lattice_OnSite
    use M_Grid

    complex(R64), intent(out), allocatable :: interactionPotential(:)
    complex(R64), intent(in), contiguous   :: src(:)
    real(R64), intent(in)                  :: time
    integer(I32), intent(in), optional     :: bt1_
    integer(I32), intent(in), optional     :: bt2_

    if (.false.) call UnusedVariables_Mark(time, bt1_, bt2_)

    if (.not. allocated(interactionPotential)) allocate (interactionPotential(Grid_nPoints))

    interactionPotential(:) = SysInteraction_Lattice_OnSite_Strength * src(:)

  end subroutine

end submodule
