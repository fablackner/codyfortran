! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for linear-grid interactions.
!>
!> @details Provides common linear-grid routines and dispatches to specific
!> interaction kernels. The source term includes quadrature weights for
!> proper numerical integration over the continuous domain.
submodule(M_SysInteraction_Linear) S_SysInteraction_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Fabricate the linear-grid interaction back-end.
  module subroutine SysInteraction_Linear_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Linear_SoftYukawa

    call Say_Fabricate("sysInteraction.linear")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_FillInteractionSrc => FillInteractionSrc
    SysInteraction_MultiplyWithInteractionPotential => MultiplyWithInteractionPotential

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysInteraction.linear.softYukawa")) then
      call SysInteraction_Linear_SoftYukawa_Fabricate

    else
      error stop "sysInteraction.linear is missing one of: softYukawa"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Build source term for linear-grid convolution.
  !>
  !> Computes the weighted density for numerical quadrature:
  !>   src(i) = conjg(orbConjg(i)) * orb(i) * w(i)
  !> where w(i) are the grid quadrature weights.
  !>
  !> @param[out] src        Allocatable complex source array (nG)
  !> @param[in]  orbConjg   Conjugated orbital (will be conjugated again → bra)
  !> @param[in]  orb        Orbital (ket)
  subroutine FillInteractionSrc(src, orbConjg, orb)
    use M_Grid
    use M_Grid_Linear

    complex(R64), intent(out), allocatable :: src(:)
    complex(R64), intent(in), contiguous :: orbConjg(:)
    complex(R64), intent(in), contiguous :: orb(:)

    integer(I32) :: nG

    nG = Grid_nPoints

    if (.not. allocated(src)) allocate (src(nG))

    src = conjg(orbConjg) * orb * Grid_Linear_weights(:)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply interaction potential to an orbital.
  !>
  !> @param[out] dOrb                   Resulting orbital dψ = V·ψ
  !> @param[in]  interactionPotential   Interaction potential V(x)
  !> @param[in]  orb                    Source orbital ψ(x)
  subroutine MultiplyWithInteractionPotential(dOrb, interactionPotential, orb)

    complex(R64), intent(out), contiguous  :: dOrb(:)
    complex(R64), intent(in), contiguous   :: interactionPotential(:)
    complex(R64), intent(in), contiguous   :: orb(:)

    dOrb(:) = interactionPotential(:) * orb(:)

  end subroutine

end submodule
