! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Dispatch submodule for MCSCF ground-solver implementation selection.
!>
!> @details
!> Reads `groundSolver.mcscf.*` and dispatches to:
!>   - `stdImpl`: Standard implementation (general grids)
!>   - `ylmOpt`: Ylm-optimized implementation (spherical grids)
!>   - `gradient`: variational MCSCF by a preconditioned orbital-gradient
!>     iteration on the MCTDHX orbital residual (any grid)
!>
!> The CI matvec is common to all of them: the CI Hamiltonian of
!> Method_Mb_OrbBased_h1/h2, which every backend fills with the current
!> orbitals before its CI diagonalization.
submodule(M_GroundSolver_Mcscf) S_GroundSolver_Mcscf

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Work arrays of GroundSolver_Mcscf_SubtractExchange: the pair potentials
  !> W(phi_p, orb) of the spin-up orbitals, their RDM-weighted sum and the
  !> product with an orbital.
  complex(R64), allocatable :: exchangePotentials(:, :)
  complex(R64), allocatable :: weightedPotential(:)
  complex(R64), allocatable :: productOrb(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Mcscf_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Method_Mb_OrbBased
    use M_GroundSolver_Mcscf_StdImpl
    use M_GroundSolver_Mcscf_YlmOpt
    use M_GroundSolver_Mcscf_Gradient

    call Say_Fabricate("groundSolver.mcscf")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    GroundSolver_Mcscf_HamiltonianAction => Method_Mb_OrbBased_ApplyCiHamiltonian

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("groundSolver.mcscf.stdImpl")) then
      call GroundSolver_Mcscf_StdImpl_Fabricate

    else if (Json_GetExistence("groundSolver.mcscf.ylmOpt")) then
      call GroundSolver_Mcscf_YlmOpt_Fabricate

    else if (Json_GetExistence("groundSolver.mcscf.gradient")) then
      call GroundSolver_Mcscf_Gradient_Fabricate

    else
      error stop "groundSolver.mcscf is missing one of: stdImpl, ylmOpt, gradient"

    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Mcscf_CiStep(coeffs, alphaCi, time)
    use M_DiagonalizerList
    use M_GroundSolver
    use M_Coeffs

    complex(R64), intent(inout), contiguous :: coeffs(:)
    real(R64), intent(in) :: alphaCi
    real(R64), intent(in) :: time

    complex(R64) :: overlap, phase

    call DiagonalizerList(1) % e % Diagonalize(time, .true.)
    call GroundSolver_NoteUnconverged(1)

    ! The eigenvector's global phase is arbitrary; align it so the mixing
    ! interpolates instead of cancelling
    overlap = dot_product(DiagonalizerList(1) % e % evecs(:, 1), coeffs)
    phase = (1.0_R64, 0.0_R64)
    if (0.0_R64 < abs(overlap)) phase = overlap / abs(overlap)

    coeffs = (1.0_R64 - alphaCi) * coeffs + alphaCi * phase * DiagonalizerList(1) % e % evecs(:, 1)
    ! The diagonalizer works on full-length vectors, whose eigenvectors lie in
    ! a declared total-M sector only up to round-off
    call Coeffs_ProjectOnSector(coeffs)
    call Coeffs_Normalize(coeffs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Mcscf_FillCorrelatedSource(weightedSrc, rdm1, orbs)
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(out), allocatable :: weightedSrc(:)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)

    complex(R64), allocatable :: src(:)
    complex(R64) :: rdm1Weight
    integer(I32) :: nOS, p, q

    nOS = size(orbs, 2)
    do p = 1, nOS
      do q = 1, nOS
        ! Restricted: the stored orbital pair carries both spins, so the
        ! density weight is the sum of the spin-up and spin-down blocks
        rdm1Weight = rdm1(p, q)
        if (Orbs_restrictedQ) rdm1Weight = rdm1Weight + rdm1(nOS + p, nOS + q)
        if (abs(rdm1Weight) < tiny(1.0_R64)) cycle
        call SysInteraction_FillInteractionSrc(src, orbs(:, p), orbs(:, q))
        if (.not. allocated(weightedSrc)) then
          allocate (weightedSrc(size(src)))
          weightedSrc = 0.0_R64
        end if
        weightedSrc = weightedSrc + rdm1Weight * src
      end do
    end do
    if (.not. allocated(weightedSrc)) error stop "GroundSolver_Mcscf_FillCorrelatedSource: the 1-RDM vanishes"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Mcscf_SubtractExchange(dOrb, orb, rdm1Up, time)
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(inout), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: orb(:)
    complex(R64), intent(in), contiguous :: rdm1Up(:, :)
    real(R64), intent(in) :: time

    complex(R64), allocatable :: src(:), interactionPotential(:)
    integer(I32) :: p, q, nUp

    nUp = size(rdm1Up, 1)

    do p = 1, nUp
      call SysInteraction_FillInteractionSrc(src, Orbs_orbs(:, p), orb(:))
      call SysInteraction_FillInteractionPotential(interactionPotential, src, time)
      if (allocated(exchangePotentials)) then
        if (size(exchangePotentials, 1) .ne. size(interactionPotential) .or. size(exchangePotentials, 2) .ne. nUp) &
          deallocate (exchangePotentials, weightedPotential)
      end if
      if (.not. allocated(exchangePotentials)) &
        allocate (exchangePotentials(size(interactionPotential), nUp), weightedPotential(size(interactionPotential)))
      exchangePotentials(:, p) = interactionPotential
    end do

    if (allocated(productOrb)) then
      if (size(productOrb) .ne. size(dOrb)) deallocate (productOrb)
    end if
    if (.not. allocated(productOrb)) allocate (productOrb(size(dOrb)))

    do q = 1, nUp
      weightedPotential = 0.0_R64
      do p = 1, nUp
        weightedPotential = weightedPotential + rdm1Up(p, q) * exchangePotentials(:, p)
      end do
      call SysInteraction_MultiplyWithInteractionPotential(productOrb, weightedPotential, Orbs_orbs(:, q))
      dOrb = dOrb - productOrb
    end do

  end subroutine

end submodule
