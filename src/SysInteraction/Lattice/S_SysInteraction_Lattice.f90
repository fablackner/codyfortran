! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for lattice-domain interactions.
!>
!> @details Provides the common lattice routines and dispatches to specific
!> interaction models (e.g., OnSite). The source term is simply the local
!> density ρᵢ = ψ*ᵢψᵢ at each site, without integration weights (discrete sum).
submodule(M_SysInteraction_Lattice) S_SysInteraction_Lattice

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Fabricate the lattice interaction back-end.
  !>
  !> Binds the common lattice routines (`FillInteractionSrc`,
  !> `MultiplyWithInteractionPotential`) and dispatches to the specific model.
  module subroutine SysInteraction_Lattice_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Lattice_OnSite

    call Say_Fabricate("sysInteraction.lattice")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_FillInteractionSrc => FillInteractionSrc
    SysInteraction_MultiplyWithInteractionPotential => MultiplyWithInteractionPotential

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysInteraction.lattice.onSite")) then
      call SysInteraction_Lattice_OnSite_Fabricate

    else
      error stop "sysInteraction.lattice is missing one of: OnSite"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Build source term for lattice interactions.
  !>
  !> Computes the local density at each lattice site:
  !>   src(i) = conjg(orbConjg(i)) * orb(i) = |ψ(i)|²
  !>
  !> @note No quadrature weights are included; the lattice is discrete.
  !>
  !> @param[out] src        Allocatable complex source array (nG)
  !> @param[in]  orbConjg   Conjugated orbital (will be conjugated again → bra)
  !> @param[in]  orb        Orbital (ket)
  subroutine FillInteractionSrc(src, orbConjg, orb)
    use M_Grid

    complex(R64), intent(out), allocatable :: src(:)
    complex(R64), intent(in), contiguous :: orbConjg(:)
    complex(R64), intent(in), contiguous :: orb(:)

    integer(I32) :: nG

    nG = Grid_nPoints

    if (.not. allocated(src)) allocate (src(nG))

    src = conjg(orbConjg) * orb
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply interaction potential to an orbital (point-wise).
  !>
  !> Performs element-wise multiplication: res(i) = V(i) * src(i)
  !>
  !> @param[out] res                    Resulting orbital
  !> @param[in]  interactionPotential   Interaction potential at each site
  !> @param[in]  src                    Source orbital
  subroutine MultiplyWithInteractionPotential(res, interactionPotential, src)

    complex(R64), intent(out), contiguous  :: res(:)
    complex(R64), intent(in), contiguous   :: interactionPotential(:)
    complex(R64), intent(in), contiguous   :: src(:)

    res(:) = interactionPotential(:) * src(:)

  end subroutine

end submodule
