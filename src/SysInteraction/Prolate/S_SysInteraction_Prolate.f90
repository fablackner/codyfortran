! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for prolate-domain interactions.
!>
!> @details Provides the common channel routines:
!>   - `FillInteractionSrc`: Product of orbitals expanded to higher mmax
!>   - `FillInteractionPotential`: Loop over m channels, solve per channel
!>   - `MultiplyWithInteractionPotential`: Channel convolution via SpatialProduct
!>
!> The channel solver (`FillInteractionPotentialChannel`) is set by the
!> specific implementation (Coulomb stdImpl).
submodule(M_SysInteraction_Prolate) S_SysInteraction_Prolate

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Fabricate the prolate interaction back-end.
  !>
  !> Reads `mmax` for the potential channels (defaults to 2x orbital mmax to
  !> capture the full azimuthal content of the density product) and dispatches
  !> to Coulomb variants.
  module subroutine SysInteraction_Prolate_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid_Prolate
    use M_SysInteraction
    use M_SysInteraction_Prolate_Coulomb

    call Say_Fabricate("sysInteraction.prolate")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Prolate_mmax = Json_Get("sysInteraction.prolate.mmax", 2 * Grid_Prolate_mmax)

    SysInteraction_FillInteractionPotential => FillInteractionPotential
    SysInteraction_FillInteractionSrc => FillInteractionSrc
    SysInteraction_MultiplyWithInteractionPotential => MultiplyWithInteractionPotential
    SysInteraction_ConjugateInteractionPotential => ConjugateInteractionPotential

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysInteraction.prolate.coulomb")) then
      call SysInteraction_Prolate_Coulomb_Fabricate

    else
      error stop "sysInteraction.prolate is missing one of: coulomb"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Compute interaction potential from source, channel by channel.
  !>
  !> Channels with negligible source (|rho_m| < 10^-14) are skipped.
  subroutine FillInteractionPotential(interactionPotential, src, time, bt1_, bt2_)
    use M_Grid
    use M_Grid_Prolate

    complex(R64), intent(out), allocatable :: interactionPotential(:)
    complex(R64), intent(in), contiguous :: src(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt1_
    integer(I32), intent(in), optional :: bt2_

    integer(I32) :: m, potSize
    complex(R64), allocatable :: srcM(:), potM(:)

    potSize = (2 * SysInteraction_Prolate_mmax + 1) * Grid_Prolate_nSpatial

    if (.not. allocated(interactionPotential)) allocate (interactionPotential(potSize))
    interactionPotential = 0.0_R64

    ! The independent channel solves may run in parallel; each channel writes
    ! a disjoint slice of the potential
    !$omp parallel default(shared) private(m, srcM, potM)
    allocate (srcM(Grid_Prolate_nSpatial), potM(Grid_Prolate_nSpatial))
    !$omp do schedule(dynamic)
    do m = -SysInteraction_Prolate_mmax, SysInteraction_Prolate_mmax
      call Grid_Prolate_GetMComponent(srcM, m, src)
      if (all(abs(srcM) < 1.0e-14_R64)) cycle
      call SysInteraction_Prolate_FillInteractionPotentialChannel(potM, srcM, m, time, bt1_, bt2_)
      call Grid_Prolate_SetMComponent(interactionPotential, m, potM)
    end do
    !$omp end do
    deallocate (srcM, potM)
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Build source term from orbital product in azimuthal channels.
  !>
  !> Computes rho = psi_i^* psi_j expanded to the potential mmax; the result
  !> includes the metric quadrature weights (Ylm convention).
  subroutine FillInteractionSrc(src, orbConjg, orb)
    use M_Grid_Prolate

    complex(R64), intent(out), allocatable :: src(:)
    complex(R64), intent(in), contiguous :: orbConjg(:)
    complex(R64), intent(in), contiguous :: orb(:)

    integer(I32) :: srcSize

    srcSize = (2 * SysInteraction_Prolate_mmax + 1) * Grid_Prolate_nSpatial

    if (.not. allocated(src)) allocate (src(srcSize))

    call Grid_Prolate_SpatialProduct(src, orbConjg, orb, SysInteraction_Prolate_mmax, &
                                     Grid_Prolate_mmax, Grid_Prolate_mmax, &
                                     conjgQ_=.true., withWeightsQ_=.true.)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply interaction potential to an orbital via channel convolution.
  subroutine MultiplyWithInteractionPotential(dOrb, interactionPotential, orb)
    use M_Grid
    use M_Grid_Prolate

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: interactionPotential(:)
    complex(R64), intent(in), contiguous :: orb(:)

    call Grid_Prolate_SpatialProduct(dOrb, interactionPotential, orb, Grid_Prolate_mmax, &
                                     SysInteraction_Prolate_mmax, Grid_Prolate_mmax)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Conjugate an interaction potential in the channel representation.
  !>
  !> Complex conjugation of a spatial function acts on the e^(i m phi)
  !> channels as f^*_m = conjg(f_(-m)).
  subroutine ConjugateInteractionPotential(potOut, potIn)
    use M_Grid_Prolate

    complex(R64), intent(out), contiguous :: potOut(:)
    complex(R64), intent(in), contiguous  :: potIn(:)

    integer(I32) :: m
    complex(R64), allocatable :: potM(:)

    potOut(:) = (0.0_R64, 0.0_R64)

    allocate (potM(Grid_Prolate_nSpatial))

    do m = -SysInteraction_Prolate_mmax, SysInteraction_Prolate_mmax
      call Grid_Prolate_GetMComponent(potM, -m, potIn)
      potM(:) = conjg(potM(:))
      call Grid_Prolate_SetMComponent(potOut, m, potM)
    end do

    deallocate (potM)

  end subroutine

end submodule
