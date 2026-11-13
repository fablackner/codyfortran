! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased) S_Method_Mb_OrbBased_Switching
  !-----------------------------------------------------------------------------
  ! Adiabatic switching (method.mb.orbBased.adiabaticSwitching.switchTime T):
  ! for -T < t < 0 the method propagates
  !   H(t) = h + (1 - lambda) U0 + lambda W,  lambda = sin^2(pi (t + T) / (2 T)),
  ! from t = 0 on the physical H = h + W. U0 is the frozen Hartree-Fock mean
  ! field of the initial state (Method_Mb_OrbBased_SetSwitchingReference): for
  ! a Hartree-Fock determinant h + U0 is its Fock operator, so at t = -T the
  ! determinant is an eigenstate of H, the closures of TD-2RDM are exact for
  ! it, and a slow switch carries it into the correlated ground state
  ! (adiabatic theorem; the gap of h + U0 protects it at the start).
  !
  ! Where the terms enter (one place each, so every method that builds its
  ! equations from these operators follows):
  !   lambda W     the source of every pair potential of S_Method_Mb_OrbBased
  !                (ApplyInteraction, FillPairPotentialsYlmNodes): h2, the
  !                correlation and Hartree-Fock operators; static orbitals
  !                scale their h2 terms by lambda (TD-2RDM-StaticOrbs);
  !   (1-lambda)U0 the one-body operator: moving orbitals on the grid
  !                (Method_Mb_OrbBased_AddSwitchingReference in
  !                ApplySingleBodyOp and the linear part), static orbitals as
  !                the matrix U0(i, j) = 2 sum_kl (h2(i,k,j,l) - h2(i,k,l,j))
  !                rdm1(l, k) of Method_Mb_OrbBased_h2 (h2 = <ij|kl>/2; also
  !                for an h2 that a program replaced) in
  !                Method_Mb_OrbBased_RefillStaticH1.
  ! Everything at t >= 0 (the setups at t = 0, static h1 and h2) is the
  ! physical Hamiltonian, so a ground-state run propagates from -T to 0 and a
  ! dynamics run continues from there. With moving orbitals the correlating
  ! orbitals start empty: a small regularizationParameter (1e-14) keeps the
  ! inverse 1-RDM from holding them back (A_Method.md, "Adiabatic Switching").
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Reference of moving orbitals (U0 on the grid) and whether it is built.
  type(T_Method_Mb_OrbBased_MeanField) :: meanField
  logical :: meanFieldQ = .false.

  !> Reference of static orbitals: the matrix of U0 (nOrbsSum x nOrbsSum).
  complex(R64), allocatable :: meanFieldMatrix(:, :)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function Method_Mb_OrbBased_SwitchFactor(time) result(res)
    use M_Utils_Constants

    real(R64), intent(in) :: time
    real(R64) :: res

    associate (switchTime => Method_Mb_OrbBased_switchTime)
      if (switchTime <= 0.0_R64 .or. time >= 0.0_R64) then
        res = 1.0_R64
      else if (time <= -switchTime) then
        res = 0.0_R64
      else
        res = sin(0.5_R64 * PI * (time + switchTime) / switchTime)**2
      end if
    end associate

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_SetSwitchingReference
    use M_Orbs
    use M_Method

    complex(R64), allocatable :: rdm1(:, :)

    if (allocated(meanFieldMatrix)) deallocate (meanFieldMatrix)
    meanFieldQ = .false.
    if (Method_Mb_OrbBased_switchTime <= 0.0_R64) return

    call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
    meanFieldQ = .not. Method_Mb_OrbBased_staticOrbsQ
    if (meanFieldQ) then
      call Method_Mb_OrbBased_FillMeanField(meanField, Orbs_orbs, rdm1, 0.0_R64)
    else
      call FillMeanFieldMatrix(meanFieldMatrix, Method_Mb_OrbBased_h2, rdm1)
    end if
    write (*, '(A,F10.4,A)') " adiabatic switching: H = h + (1 - lambda) U0 + lambda W over ", &
      Method_Mb_OrbBased_switchTime, " a.u. before t = 0, U0 the frozen mean field of the state"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_AddSwitchingReference(dOrbs, orbs, time)
    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    real(R64), intent(in) :: time

    complex(R64), allocatable :: weightedOrb(:), scratch(:)
    real(R64) :: weight
    integer(I32) :: iOrb

    if (.not. meanFieldQ) return
    weight = 1.0_R64 - Method_Mb_OrbBased_SwitchFactor(time)
    if (weight <= 0.0_R64) return

    ! serial over the orbitals, parallel inside
    allocate (weightedOrb(size(orbs, 1)), scratch(size(orbs, 1)))
    do iOrb = 1, size(orbs, 2)
      weightedOrb(:) = weight * orbs(:, iOrb)
      call Method_Mb_OrbBased_AddMeanField(dOrbs(:, iOrb), scratch, meanField, weightedOrb, &
                                           Method_Mb_OrbBased_bodyTypeOfOrb(iOrb))
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_AddSwitchingReferenceMatrix(h1, time)
    complex(R64), intent(inout), contiguous :: h1(:, :)
    real(R64), intent(in) :: time

    real(R64) :: weight

    if (.not. allocated(meanFieldMatrix)) return
    weight = 1.0_R64 - Method_Mb_OrbBased_SwitchFactor(time)
    if (weight <= 0.0_R64) return

    h1(:, :) = h1(:, :) + weight * meanFieldMatrix(:, :)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Mean-field matrix of h2 and rdm1 (fermions, h2 = <ij|kl>/2, zero for
  !> pairs of different body type, rdm1(i, j) = <a+_j a_i>): the derivative of
  !> the energy of the determinant, U(i, j) = 2 sum_kl (h2(i,k,j,l) -
  !> h2(i,k,l,j)) rdm1(l, k).
  subroutine FillMeanFieldMatrix(u, h2, rdm1)
    complex(R64), intent(out), allocatable :: u(:, :)
    complex(R64), intent(in), contiguous :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous :: rdm1(:, :)

    integer(I32) :: n, i, j, k, l

    n = size(rdm1, 1)
    allocate (u(n, n))
    !$omp parallel do default(shared) private(i, j, k, l) collapse(2)
    do j = 1, n
      do i = 1, n
        u(i, j) = (0.0_R64, 0.0_R64)
        do l = 1, n
          do k = 1, n
            u(i, j) = u(i, j) + 2.0_R64 * (h2(i, k, j, l) - h2(i, k, l, j)) * rdm1(l, k)
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine

end submodule
