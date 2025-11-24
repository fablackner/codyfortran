! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_GroundSolver_Tdhx_StdImpl) S_GroundSolver_Tdhx_StdImpl

  implicit none

  complex(R64), allocatable :: dOrbTmp(:)
  complex(R64), allocatable :: hartreePotential(:)
  complex(R64), allocatable :: externalPotential(:)
  complex(R64), allocatable :: interactionPotential(:)
  complex(R64), allocatable :: src(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Tdhx_StdImpl_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_GroundSolver
    use M_GroundSolver_Tdhx

    call Say_Fabricate("groundSolver.tdhx.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    GroundSolver_Setup => Setup
    GroundSolver_Approach => Approach
    GroundSolver_Tdhx_HartreeFockAction => HartreeFockAction

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm
    use M_SysInteraction_Ylm

    integer(I32) :: lmaxPot, potSize

    call Say_Setup("groundSolver.tdhx.stdImpl")

    lmaxPot = SysInteraction_Ylm_lmax
    potSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    allocate (hartreePotential(potSize))
    allocate (dOrbTmp(Grid_nPoints))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine HartreeFockAction(dOrb, orb, time)
    use M_SysKinetic
    use M_SysPotential
    use M_SysInteraction
    use M_Method_Mb_OrbBased
    use M_Orbs

    complex(R64), intent(out), contiguous, target :: dOrb(:)
    complex(R64), intent(in), contiguous, target :: orb(:)
    real(R64), intent(in) :: time

    integer(I32) :: j

    dOrb = 0.0_R64

    ! Kinetic part
    call SysKinetic_MultiplyWithKineticOp(dOrbTmp, orb, time)
    dOrb = dOrb + dOrbTmp

    ! External potential part
    call SysPotential_MultiplyWithExternalPotential(dOrbTmp, externalPotential, orb)
    dOrb = dOrb + dOrbTmp

    ! Interaction part (Hartree-Fock)
    call SysInteraction_MultiplyWithInteractionPotential(dOrbTmp, hartreePotential, orb)
    dOrb = dOrb + dOrbTmp

    do j = 1, Orbs_nOrbsInState / 2
      call SysInteraction_FillInteractionSrc(src, Orbs_orbs(:, j), orb(:))
      call SysInteraction_FillInteractionPotential(interactionPotential, src, time)
      call SysInteraction_MultiplyWithInteractionPotential(dOrbTmp, interactionPotential, Orbs_orbs(:, j))
      dOrb = dOrb - dOrbTmp
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Approach(state, alpha, time)
    use M_Grid
    use M_Grid_Ylm
    use M_DiagonalizerList
    use M_Orbs
    use M_SysPotential
    use M_SysInteraction
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: alpha
    real(R64), intent(in) :: time

    complex(R64), contiguous, pointer :: orbs(:, :)
    integer(I32) :: nG, nOS, j

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    orbs(1:nG, 1:nOS) => state(1:)

    call SysPotential_FillExternalPotential(externalPotential, time)

    hartreePotential = 0.0_R64
    do j = 1, nOS / 2
      ! Direct (Hartree) term assume body-type independence
      call SysInteraction_FillInteractionSrc(src, orbs(:, j), orbs(:, j))
      call SysInteraction_FillInteractionPotential(interactionPotential, 2.0_R64 * src, time)
      hartreePotential = hartreePotential + interactionPotential
    end do

    ! Diagonalize
    call DiagonalizerList(1) % e % Diagonalize(time, .true.)

    ! Mix the new orbitals with the old ones
    do j = 1, nOS / 2
      orbs(:, j) = (1.0_R64 - alpha) * orbs(:, j) + alpha * DiagonalizerList(1) % e % evecs(:, j)
      orbs(:, nOS / 2 + j) = orbs(:, j)
    end do

    ! Orthonormalize the orbitals
    call Orbs_Orthonormalize(orbs)

  end subroutine

end submodule
