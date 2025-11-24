! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_GroundSolver_Tdhx_YlmOpt) S_GroundSolver_Tdhx_YlmOpt

  implicit none

  complex(R64), allocatable :: dOrbLmTmp(:)
  complex(R64), allocatable :: potLm(:)
  complex(R64), allocatable :: potLmTmp(:)
  complex(R64), allocatable :: src(:)
  complex(R64), allocatable :: srcTmp(:)
  complex(R64), allocatable :: srcLm(:)
  complex(R64), allocatable :: interactionPotential(:)
  complex(R64), allocatable :: orb(:)
  complex(R64), allocatable :: dOrbTmp(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Tdhx_YlmOpt_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_GroundSolver
    use M_GroundSolver_Tdhx

    call Say_Fabricate("groundSolver.tdhx.ylmOpt")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    GroundSolver_Setup => Setup
    GroundSolver_Approach => Approach
    GroundSolver_Tdhx_YlmOpt_HartreeFockAction => HartreeFockAction

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm
    use M_SysInteraction_Ylm

    integer(I32) :: lmaxPot, potSize

    call Say_Setup("groundSolver.tdhx.ylmOpt")

    lmaxPot = SysInteraction_Ylm_lmax
    potSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    allocate (src(potSize))
    allocate (srcTmp(potSize))
    allocate (interactionPotential(potSize))
    allocate (dOrbLmTmp(Grid_Ylm_nRadial))
    allocate (potLm(Grid_Ylm_nRadial))
    allocate (potLmTmp(Grid_Ylm_nRadial))
    allocate (srcLm(Grid_Ylm_nRadial))
    allocate (orb(Grid_nPoints))
    allocate (dOrbTmp(Grid_nPoints))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine HartreeFockAction(dOrbLm, orbLm, l, time)
    use M_Utils_SphericalHarmonics
    use M_SysKinetic_Ylm
    use M_SysPotential_Ylm
    use M_SysInteraction
    use M_Method_Mb_OrbBased
    use M_Orbs
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous, target :: dOrbLm(:)
    complex(R64), intent(in), contiguous, target :: orbLm(:)
    integer(I32), intent(in) :: l
    real(R64), intent(in) :: time

    real(R64) :: gVal
    integer(I32) :: j

    dOrbLm = 0.0_R64

    gVal = SphericalHarmonics_GauntCoefficient( &
           l1=0, m1=0, &
           l2=l, m2=0, &
           l3=l, m3=0)

    ! Kinetic part
    call SysKinetic_Ylm_MultiplyWithRadialKineticOp(dOrbLmTmp, orbLm, l, 0, time)
    dOrbLm = dOrbLm + dOrbLmTmp

    ! Potential part
    dOrbLm = dOrbLm + gVal * potLm * orbLm

    orb = 0.0_R64
    call Grid_Ylm_SetLmComponent(orb, l, 0, orbLm)
    do j = 1, Orbs_nOrbsInState / 2
      call SysInteraction_FillInteractionSrc(src, Orbs_orbs(:, j), orb(:))
      call SysInteraction_FillInteractionPotential(interactionPotential, src, time)
      call SysInteraction_MultiplyWithInteractionPotential(dOrbTmp, interactionPotential, Orbs_orbs(:, j))
      call Grid_Ylm_GetLmComponent(dOrbLmTmp, l, 0, dOrbTmp)
      dOrbLm = dOrbLm - dOrbLmTmp
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Approach(state, alpha, time)
    use M_Grid
    use M_Grid_Ylm
    use M_DiagonalizerList
    use M_Orbs
    use M_SysInteraction
    use M_SysInteraction_Ylm
    use M_SysPotential_Ylm
    use M_OrbsInit_Ylm_HydrogenLike

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: alpha
    real(R64), intent(in) :: time

    complex(R64), contiguous, pointer :: orbs(:, :)
    integer(I32) :: nG, nOS, j
    integer(I32) :: n, l, m, nr

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    orbs(1:nG, 1:nOS) => state(1:)

    call SysPotential_Ylm_FillExternalPotentialRadial(potLm, 0, 0, time)

    src = 0.0_R64
    do j = 1, nOS / 2
      call SysInteraction_FillInteractionSrc(srcTmp, orbs(:, j), orbs(:, j))
      src = src + 2.0_R64 * srcTmp
    end do
    call Grid_Ylm_GetLmComponent(srcLm, 0, 0, src)
    call SysInteraction_Ylm_FillInteractionPotentialRadial(potLmTmp, srcLm, 0, 0, time)
    potLm = potLm + potLmTmp

    ! Diagonalize
    do l = 0, Grid_Ylm_lmax
      call DiagonalizerList(l + 1) % e % Diagonalize(time, .true.)
    end do

    do j = 1, nOS / 2
      n = OrbsInit_Ylm_HydrogenLike_n(j)
      l = OrbsInit_Ylm_HydrogenLike_l(j)
      m = OrbsInit_Ylm_HydrogenLike_m(j)
      nr = n - l - 1

      dOrbTmp = 0.0_R64
      call Grid_Ylm_SetLmComponent(dOrbTmp(:), l, m, DiagonalizerList(l + 1) % e % evecs(:, nr + 1))
      orbs(:, j) = (1.0_R64 - alpha) * orbs(:, j) + alpha * dOrbTmp(:)
      orbs(:, nOS / 2 + j) = orbs(:, j)
    end do

    ! Orthonormalize the mixed orbitals
    call Orbs_Orthonormalize(orbs)

  end subroutine

end submodule
