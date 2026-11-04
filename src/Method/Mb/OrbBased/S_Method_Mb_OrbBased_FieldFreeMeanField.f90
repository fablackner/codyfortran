! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased) S_Method_Mb_OrbBased_FieldFreeMeanField
  !-----------------------------------------------------------------------------
  ! Frozen mean-field reference of the exact field-free part
  ! (method.mb.orbBased.exactFieldFreeMeanFieldQ, MCTDHX with the Cayley
  ! exact part of S_Method_Mb_OrbBased_FieldFreeCayley): a static one-body
  ! operator U = J0 - K0 of the state at the first propagation
  ! (Method_Mb_OrbBased_BuildFieldFreeMeanField, from its 1-RDM): J0 the
  ! Hartree potential of its density (SysInteraction_FillDirectPotential:
  ! one potential solve), K0 its exchange operator compressed with ACE
  ! (SysInteraction_FillExchangeAce, n(n+1)/2 pair potentials once),
  ! K0 v ~ sum_k zeta_k <zeta_k|v>, exact on the span of the orbitals and
  ! Hermitian, so A = h0 + U stays Hermitian in the grid metric and the
  ! Cayley step unitary. U v is SysInteraction_AddMeanField with J0 and the
  ! ACE vectors of the body type of v (the Fock action of the SCF with aceQ
  ! uses the same helper).
  !
  ! Two users, each with its own entry point: the Cayley step propagates A
  ! and applies U per vector in its Arnoldi operator and true residual
  ! (Method_Mb_OrbBased_ApplyFieldFreeMeanFieldOnOrb: one multiplication,
  ! n inner products and one fused pass per application; the preconditioner
  ! is unchanged), and MCTDHX subtracts U in its nonlinear part
  ! (Method_Mb_OrbBased_ApplyFieldFreeMeanField), so the equations of motion
  ! are unchanged. Why: in a heavy core (the O 1s of H2O) h0 phi and the mean
  ! field W phi are both large and cancel only in their sum, so the split
  ! h0 | W has a large error; the split (h0 + U) | (W - U) has a remainder of
  ! the order of the deviation from the reference state (O(kappa) after a
  ! kick).
  !
  ! Lifetime: every setup of the Cayley part ends with
  ! Method_Mb_OrbBased_SetupFieldFreeMeanField, which releases U; the next
  ! propagation builds it again unless a program does.
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Whether Method_Mb_OrbBased_SetupFieldFreeMeanField has run.
  logical :: setupQ = .false.

  !> Frozen mean-field reference U = J0 - K0 (exactFieldFreeMeanFieldQ,
  !> Method_Mb_OrbBased_BuildFieldFreeMeanField): whether it is built since
  !> the last setup, the Hartree potential J0 of the reference density (in
  !> the representation of the interaction backend) and the ACE vectors of
  !> K0, one per stored orbital; the columns meanFieldZetaStart(bt) ..
  !> meanFieldZetaEnd(bt) act on the orbitals of body type bt.
  logical :: meanFieldBuiltQ = .false.
  complex(R64), allocatable :: meanFieldHartree(:)
  complex(R64), allocatable :: meanFieldZeta(:, :)
  integer(I32), allocatable :: meanFieldZetaStart(:)
  integer(I32), allocatable :: meanFieldZetaEnd(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_SetupFieldFreeMeanField
    !---------------------------------------------------------------------------
    ! Releases the reference of an earlier setup (the next propagation builds
    ! it for this one unless a program does); with exactFieldFreeMeanFieldQ
    ! the interaction must be body-type independent. Called at the end of
    ! Method_Mb_OrbBased_SetupFieldFreeCayley, so its line follows the GMRES
    ! line there.
    !---------------------------------------------------------------------------
    use M_SysInteraction

    call ReleaseMeanField
    if (Method_Mb_OrbBased_exactFieldFreeMeanFieldQ) then
      if (.not. SysInteraction_bodyTypeIndependentQ) &
        error stop "method.mb.orbBased.exactFieldFreeMeanFieldQ requires a body-type-independent interaction &
                   &(sysInteraction)"
    end if
    setupQ = .true.

    if (Method_Mb_OrbBased_exactFieldFreeMeanFieldQ) &
      write (*, '(A)') " exact field-free part: h0 + U with the frozen mean-field reference U = J0 - K0 &
                       &(exactFieldFreeMeanFieldQ), built from the state at the first propagation unless built before"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_BuildFieldFreeMeanField(orbs, rdm1, time)
    !---------------------------------------------------------------------------
    ! U = J0 - K0 from the stored orbitals phi and the 1-RDM of the state
    ! (rdm1(i, j) = <a+_j a_i>), the factors of ApplyHartreeFockOp:
    !   J0 = V[rho], rho = sum_kl rhoS(k, l) conjg(phi_l) phi_k with the
    !        spatial 1-RDM rhoS (Method_Mb_OrbBased_FillSpatialRdm1: the sum
    !        of both spin blocks for restricted storage, i.e. 2 sum |phi|^2
    !        for a closed shell), one solve of the total source;
    !   K0 within every body type bt (exchange only between equal spins):
    !        gamma_bt(r, r') = sum_{j,l in bt} rdm1(j, l) phi_j(r) conjg(phi_l(r')),
    !        SysInteraction_FillExchangeAce on the stored orbitals of bt with
    !        the weighted orbitals phi rdm1_bt (restricted storage: one block,
    !        all stored orbitals with the spin-up block of rdm1, shared by
    !        both spins as the orbitals are).
    ! The potentials take the body-type-independent arguments (-1, -1) of
    ! ApplyInteraction (checked at setup).
    !---------------------------------------------------------------------------
    use M_Grid
    use M_Orbs
    use M_Method_Mb
    use M_SysInteraction

    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    real(R64), intent(in) :: time

    complex(R64), allocatable :: rdmSpatial(:, :), weightedOrbs(:, :)
    real(R64) :: seconds
    integer(I64) :: count0, count1, countRate
    integer(I32) :: nOS, nPairSolves, bt, iStart, iEnd, n
    logical :: positiveDefiniteQ

    if (.not. setupQ) &
      error stop "method.mb.orbBased.exactFieldFreeMeanFieldQ: call Method_Setup before building the reference"
    nOS = Orbs_nOrbsInState
    if (size(orbs, 1) .ne. Grid_nPoints .or. size(orbs, 2) .ne. nOS) &
      error stop "Method_Mb_OrbBased_BuildFieldFreeMeanField: orbs must be the stored orbitals (nPoints x nOrbsInState)"
    if (size(rdm1, 1) .ne. Method_Mb_OrbBased_nOrbsSum .or. size(rdm1, 2) .ne. Method_Mb_OrbBased_nOrbsSum) &
      error stop "Method_Mb_OrbBased_BuildFieldFreeMeanField: rdm1 must be nOrbsSum x nOrbsSum"

    call system_clock(count0, countRate)
    call ReleaseMeanField

    ! J0: one potential of the total source sum_l conjg(phi_l) (phi rhoS)_l
    call Method_Mb_OrbBased_FillSpatialRdm1(rdmSpatial, rdm1)
    weightedOrbs = matmul(orbs, rdmSpatial)
    call SysInteraction_FillDirectPotential(meanFieldHartree, orbs, weightedOrbs, time, bt1_=-1, bt2_=-1)

    ! K0: ACE vectors per body type of the stored orbitals
    allocate (meanFieldZeta(Grid_nPoints, nOS))
    allocate (meanFieldZetaStart(Method_Mb_nBodyTypes), meanFieldZetaEnd(Method_Mb_nBodyTypes))
    nPairSolves = 0
    do bt = 1, Method_Mb_nBodyTypes
      if (Orbs_restrictedQ) then
        iStart = 1
        iEnd = nOS
      else
        iStart = Method_Mb_OrbBased_nOrbsStart(bt)
        iEnd = Method_Mb_OrbBased_nOrbsEnd(bt)
      end if
      meanFieldZetaStart(bt) = iStart
      meanFieldZetaEnd(bt) = iEnd

      ! Restricted: one block (spin-up rdm1) for the shared orbitals
      if (Orbs_restrictedQ .and. bt > 1) cycle
      n = iEnd - iStart + 1
      if (n < 1) cycle
      weightedOrbs = matmul(orbs(:, iStart:iEnd), rdm1(iStart:iEnd, iStart:iEnd))
      call SysInteraction_FillExchangeAce(meanFieldZeta(:, iStart:iEnd), positiveDefiniteQ, orbs(:, iStart:iEnd), &
                                          weightedOrbs, time, bt1_=-1, bt2_=-1)
      if (.not. positiveDefiniteQ) &
        write (*, '(A,I0,A)') " warning: exactFieldFreeMeanFieldQ: the exchange matrix <phi_i|K0 phi_j> of body type ", bt, &
        " is not positive definite, its ACE uses the positive eigen-subspace"
      if (SysInteraction_conjSymmetricQ) then
        nPairSolves = nPairSolves + n * (n + 1) / 2
      else
        nPairSolves = nPairSolves + n * n
      end if
    end do

    meanFieldBuiltQ = .true.

    call system_clock(count1)
    seconds = real(count1 - count0, kind=R64) / real(countRate, kind=R64)
    write (*, '(A,F10.4,A,F14.10,A,I0,A,F8.2,A)') " exact field-free part: frozen mean-field reference U = J0 - K0 built at t =", &
      time, " from the 1-RDM of the state (N =", real(SumDiagonal(rdmSpatial), kind=R64), "): 1 + ", nPairSolves, &
      " potential solves,", seconds, " s"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_ApplyFieldFreeMeanField(dOrbs, orbs, h1_)
    !---------------------------------------------------------------------------
    ! dOrbs(:, i) += U orbs(:, i) (body type of orbital i), and optionally
    ! h1_(i, j) = <orbs_i|U orbs_j> for equal body types (restricted: the
    ! spatial block in both spin blocks), as Method_Mb_OrbBased_FillH1
    ! (Method_Mb_OrbBased_FillOneBodyMatrix). Serial over the orbitals,
    ! parallel inside (call it outside parallel regions for speed; it is
    ! reentrant).
    !---------------------------------------------------------------------------
    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(out), allocatable, optional :: h1_(:, :)

    complex(R64), allocatable :: uOrbs(:, :), scratch(:)
    integer(I32) :: nOS, iOrb

    if (.not. meanFieldBuiltQ) &
      error stop "Method_Mb_OrbBased_ApplyFieldFreeMeanField: the frozen mean-field reference is not built &
                 &(Method_Mb_OrbBased_BuildFieldFreeMeanField)"

    nOS = size(orbs, 2)
    allocate (uOrbs(size(orbs, 1), nOS), scratch(size(orbs, 1)))
    do iOrb = 1, nOS
      uOrbs(:, iOrb) = (0.0_R64, 0.0_R64)
      call AddMeanField(uOrbs(:, iOrb), scratch, orbs(:, iOrb), Method_Mb_OrbBased_bodyTypeOfOrb(iOrb))
      call Axpy(dOrbs(:, iOrb), (1.0_R64, 0.0_R64), uOrbs(:, iOrb))
    end do

    if (present(h1_)) call Method_Mb_OrbBased_FillOneBodyMatrix(h1_, orbs, uOrbs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_ApplyFieldFreeMeanFieldOnOrb(dOrb, scratch, orb, bt)
    !---------------------------------------------------------------------------
    ! dOrb += U orb of body type bt (AddMeanField) once the reference is built
    ! and exactFieldFreeMeanFieldQ is set (the MCTDHX nonlinear part subtracts
    ! U under the same flag); a no-op otherwise.
    !---------------------------------------------------------------------------
    complex(R64), intent(inout), contiguous :: dOrb(:)
    complex(R64), intent(inout), contiguous :: scratch(:)
    complex(R64), intent(in), contiguous :: orb(:)
    integer(I32), intent(in) :: bt

    if (meanFieldBuiltQ .and. Method_Mb_OrbBased_exactFieldFreeMeanFieldQ) call AddMeanField(dOrb, scratch, orb, bt)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function Method_Mb_OrbBased_FieldFreeMeanFieldBuiltQ() result(res)
    logical :: res

    res = meanFieldBuiltQ

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> w = w + U z = w + J0 z - sum_k zeta_k <zeta_k|z> with the ACE vectors of
  !> body type bt (SysInteraction_AddMeanField; scratch: work vector for
  !> J0 z). Reentrant; its passes are parallel when called outside a parallel
  !> region.
  subroutine AddMeanField(w, scratch, z, bt)
    use M_SysInteraction

    complex(R64), intent(inout), contiguous :: w(:)
    complex(R64), intent(inout), contiguous :: scratch(:)
    complex(R64), intent(in), contiguous :: z(:)
    integer(I32), intent(in) :: bt

    call SysInteraction_AddMeanField(w, scratch, z, meanFieldHartree, &
                                     meanFieldZeta(:, meanFieldZetaStart(bt):meanFieldZetaEnd(bt)))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Releases the frozen mean-field reference.
  subroutine ReleaseMeanField

    meanFieldBuiltQ = .false.
    if (allocated(meanFieldHartree)) deallocate (meanFieldHartree)
    if (allocated(meanFieldZeta)) deallocate (meanFieldZeta)
    if (allocated(meanFieldZetaStart)) deallocate (meanFieldZetaStart)
    if (allocated(meanFieldZetaEnd)) deallocate (meanFieldZetaEnd)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Sum of the diagonal of a square matrix.
  pure function SumDiagonal(a) result(res)
    complex(R64), intent(in) :: a(:, :)
    complex(R64) :: res

    integer(I32) :: i

    res = (0.0_R64, 0.0_R64)
    do i = 1, min(size(a, 1), size(a, 2))
      res = res + a(i, i)
    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> y = y + a x.
  !> Private twin of the Axpy of M_Utils_Gmres, kept for bit-identity (the former (1,0)*u arithmetic).
  subroutine Axpy(y, a, x)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(inout), contiguous :: y(:)
    complex(R64), intent(in) :: a
    complex(R64), intent(in), contiguous :: x(:)

    integer(I32) :: i

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(y)
      y(i) = y(i) + a * x(i)
    end do
    !$omp end parallel do

  end subroutine

end submodule
