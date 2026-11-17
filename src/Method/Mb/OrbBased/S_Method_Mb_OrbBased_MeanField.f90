! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased) S_Method_Mb_OrbBased_MeanField
  !-----------------------------------------------------------------------------
  ! Frozen Hartree-Fock mean field U = J0 - K0 of a 1-RDM on the grid
  ! (T_Method_Mb_OrbBased_MeanField): J0 the Hartree potential of its density
  ! (SysInteraction_FillDirectPotential: one potential solve), K0 its exchange
  ! operator compressed with ACE (SysInteraction_FillExchangeAce, n(n+1)/2
  ! pair potentials once), K0 v ~ sum_k zeta_k <zeta_k|v>, exact on the span
  ! of the orbitals and Hermitian. U v is SysInteraction_AddMeanField with J0
  ! and the ACE vectors of the body type of v (the Fock action of the SCF with
  ! aceQ uses the same helper).
  !
  ! Two users: the frozen reference of the Cayley exact part
  ! (S_FieldFree_Cayley_MeanField) and the reference of the adiabatic
  ! switching (S_Method_Mb_OrbBased_Switching).
  !-----------------------------------------------------------------------------

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_FillMeanField(mf, orbs, rdm1, time, nPairSolves_, exchangeOnlyQ_)
    !---------------------------------------------------------------------------
    ! U = J0 - K0 from the stored orbitals phi and the 1-RDM of the state
    ! (rdm1(i, j) = <a+_j a_i>), the factors of ApplyHartreeFockOp:
    !   J0 = V[rho], rho = sum_kl rhoS(k, l) conjg(phi_l) phi_k with the
    !        spatial 1-RDM rhoS (Method_Mb_OrbBased_FillSpatialRdm1: the sum
    !        of both spin blocks for restricted storage, i.e. 2 sum |phi|^2
    !        for a closed shell), one solve of the total source;
    !   K0 by its ACE vectors (FillExchangeAce).
    ! The potentials take the body-type-independent arguments (-1, -1) of
    ! ApplyInteraction. With exchangeOnlyQ_ only the ACE vectors are built
    ! (mf % hartree stays unallocated).
    !---------------------------------------------------------------------------
    use M_SysInteraction

    type(T_Method_Mb_OrbBased_MeanField), intent(out) :: mf
    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    real(R64), intent(in) :: time
    integer(I32), intent(out), optional :: nPairSolves_
    logical, intent(in), optional :: exchangeOnlyQ_

    complex(R64), allocatable :: rdmSpatial(:, :)
    integer(I32) :: nPairSolves
    logical :: exchangeOnlyQ

    exchangeOnlyQ = .false.
    if (present(exchangeOnlyQ_)) exchangeOnlyQ = exchangeOnlyQ_

    if (.not. SysInteraction_bodyTypeIndependentQ) &
      error stop "Method_Mb_OrbBased_FillMeanField requires a body-type-independent interaction (sysInteraction)"

    if (.not. exchangeOnlyQ) then
      call Method_Mb_OrbBased_FillSpatialRdm1(rdmSpatial, rdm1)
      call SysInteraction_FillDirectPotential(mf % hartree, orbs, matmul(orbs, rdmSpatial), time, bt1_=-1, bt2_=-1)
    end if
    call FillExchangeAce(mf, nPairSolves, orbs, rdm1, time)
    if (present(nPairSolves_)) nPairSolves_ = nPairSolves

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_AddMeanField(w, scratch, mf, z, bt)
    use M_SysInteraction

    complex(R64), intent(inout), contiguous :: w(:)
    complex(R64), intent(inout), contiguous :: scratch(:)
    type(T_Method_Mb_OrbBased_MeanField), intent(in) :: mf
    complex(R64), intent(in), contiguous :: z(:)
    integer(I32), intent(in) :: bt

    call SysInteraction_AddMeanField(w, scratch, z, mf % hartree, mf % zeta(:, mf % zetaStart(bt):mf % zetaEnd(bt)))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> The ACE vectors of K0 within every body type bt (exchange only between
  !> equal spins): gamma_bt(r, r') = sum_{j,l in bt} rdm1(j, l) phi_j(r)
  !> conjg(phi_l(r')), SysInteraction_FillExchangeAce on the stored orbitals
  !> of bt with the weighted orbitals phi rdm1_bt (restricted storage: one
  !> block, all stored orbitals with the spin-up block of rdm1, shared by both
  !> spins as the orbitals are); the columns zetaStart(bt) .. zetaEnd(bt) act
  !> on the orbitals of body type bt. nPairSolves: the pair potentials solved.
  subroutine FillExchangeAce(mf, nPairSolves, orbs, rdm1, time)
    use M_Grid
    use M_Orbs
    use M_Method_Mb
    use M_SysInteraction

    type(T_Method_Mb_OrbBased_MeanField), intent(inout) :: mf
    integer(I32), intent(out) :: nPairSolves
    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    real(R64), intent(in) :: time

    integer(I32) :: nOS, bt, iStart, iEnd, n
    logical :: positiveDefiniteQ

    nOS = Orbs_nOrbsInState
    allocate (mf % zeta(Grid_nPoints, nOS))
    allocate (mf % zetaStart(Method_Mb_nBodyTypes), mf % zetaEnd(Method_Mb_nBodyTypes))
    nPairSolves = 0
    do bt = 1, Method_Mb_nBodyTypes
      if (Orbs_restrictedQ) then
        iStart = 1
        iEnd = nOS
      else
        iStart = Method_Mb_OrbBased_nOrbsStart(bt)
        iEnd = Method_Mb_OrbBased_nOrbsEnd(bt)
      end if
      mf % zetaStart(bt) = iStart
      mf % zetaEnd(bt) = iEnd

      ! Restricted: one block (spin-up rdm1) for the shared orbitals
      if (Orbs_restrictedQ .and. bt > 1) cycle
      n = iEnd - iStart + 1
      if (n < 1) cycle
      call SysInteraction_FillExchangeAce(mf % zeta(:, iStart:iEnd), positiveDefiniteQ, orbs(:, iStart:iEnd), &
                                          matmul(orbs(:, iStart:iEnd), rdm1(iStart:iEnd, iStart:iEnd)), time, &
                                          bt1_=-1, bt2_=-1)
      if (.not. positiveDefiniteQ) &
        write (*, '(A,I0,A)') " warning: frozen mean field: the exchange matrix <phi_i|K0 phi_j> of body type ", bt, &
        " is not positive definite, its ACE uses the positive eigen-subspace"
      if (SysInteraction_conjSymmetricQ) then
        nPairSolves = nPairSolves + n * (n + 1) / 2
      else
        nPairSolves = nPairSolves + n * n
      end if
    end do

  end subroutine

end submodule
