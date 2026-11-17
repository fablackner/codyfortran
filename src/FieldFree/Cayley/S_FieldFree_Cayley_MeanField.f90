! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_FieldFree_Cayley:S_FieldFree_Cayley) S_FieldFree_Cayley_MeanField
!-----------------------------------------------------------------------------
! Frozen mean-field reference of the Cayley exact part
! (fieldFree.cayley.meanFieldQ): a static one-body operator U = J0 - K0 of
! the state at the first propagation (FieldFree_BuildMeanField, from its
! 1-RDM), the mean field of S_Method_Mb_OrbBased_MeanField
! (Method_Mb_OrbBased_FillMeanField): Hermitian, so A = h0 + U stays
! Hermitian in the grid metric and the Cayley step unitary.
!
! A child of S_FieldFree_Cayley (nested submodule): the reference is an
! add-on of the Cayley step, which declares the procedures implemented here
! (SetupMeanField, ApplyMeanFieldOnOrb and the targets BuildMeanField and
! ApplyMeanField of FieldFree_BuildMeanField and FieldFree_ApplyMeanField).
!
! Two users, each with its own entry point: the Cayley step propagates A
! and applies U per vector in its Arnoldi operator and true residual
! (ApplyMeanFieldOnOrb: one multiplication, n inner products and one fused
! pass per application; the preconditioner is unchanged), and the method
! subtracts U in its nonlinear part (FieldFree_ApplyMeanField), so the
! equations of motion are unchanged. Why: in a heavy core (the O 1s of H2O) h0 phi and the mean
! field W phi are both large and cancel only in their sum, so the split
! h0 | W has a large error; the split (h0 + U) | (W - U) has a remainder of
! the order of the deviation from the reference state (O(kappa) after a
! kick).
!
! Lifetime: every setup of the Cayley part ends with SetupMeanField, which
! releases U; the next propagation builds it again unless a program does.
!
! Kohn-Sham reference (Method_kohnShamQ, method.mb.orbBased.tdks):
! U = v_KS,s[rho0] - a K0 with the local Kohn-Sham potential of the state on
! the xc points of SysXc (one column per spin channel,
! Method_Mb_OrbBased_Tdks_FillKohnShamPotential), applied by
! SysXc_MultiplyWithPotential, and the ACE vectors of K0 scaled by
! sqrt(a) (a = SysXc_exactExchangeFraction; none for a = 0), so that
! sum_k zeta_k <zeta_k|v> = a K0 v on the span of the orbitals. The TDKS
! nonlinear part subtracts U: at the Kohn-Sham ground state it vanishes.
!-----------------------------------------------------------------------------

use M_FieldFree, only: FieldFree_meanFieldQ, FieldFree_meanFieldBuiltQ
use M_Method_Mb_OrbBased

implicit none

!=============================================================================
! local data
!=============================================================================

!> Whether SetupMeanField has run.
logical :: referenceSetupQ = .false.

!> Frozen mean-field reference U (FieldFree_BuildMeanField; built or not:
!> FieldFree_meanFieldBuiltQ): J0 and the ACE vectors of K0 (Hartree-Fock, or
!> only the ACE vectors of a K0 for Kohn-Sham).
type(T_Method_Mb_OrbBased_MeanField) :: ref

!> Kohn-Sham reference: the local potential v_KS,s[rho0] on the xc points
!> (SysXc_nPoints, nSpin), the column of every body type, and whether ACE
!> vectors (a K0, a > 0) exist.
real(R64), allocatable :: meanFieldKohnSham(:, :)
integer(I32), allocatable :: meanFieldSpinOfBodyType(:)
logical :: meanFieldExchangeQ = .true.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module subroutine SetupMeanField
  !---------------------------------------------------------------------------
  ! Releases the reference of an earlier setup (the next propagation builds
  ! it for this one unless a program does); with FieldFree_meanFieldQ the
  ! bodies must be fermions and the interaction body-type independent. Called
  ! at the end of the Cayley setup, so its line follows the GMRES line there.
  !---------------------------------------------------------------------------
  use M_Utils_Json
  use M_SysInteraction
  use M_Method_Mb

  call ReleaseMeanField
  if (FieldFree_meanFieldQ) then
    if (any(Method_Mb_bodyStatistics .ne. 'f')) &
      error stop "fieldFree.cayley.meanFieldQ requires fermions (the reference is the Hartree-Fock or Kohn-Sham &
                 &mean field of the 1-RDM)"
    if (Json_GetExistence("sysInteraction")) then
      if (.not. SysInteraction_bodyTypeIndependentQ) &
        error stop "fieldFree.cayley.meanFieldQ requires a body-type-independent interaction (sysInteraction)"
    end if
  end if
  referenceSetupQ = .true.

  if (FieldFree_meanFieldQ) &
    write (*, '(3A)') " exact field-free part: h0 + U with the ", ReferenceName(), &
    " (fieldFree.cayley.meanFieldQ), built from the state at the first propagation unless built before"

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module subroutine BuildMeanField(orbs, rdm1, time)
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
  use M_SysXc
  use M_Method, only: Method_kohnShamQ
  use M_Method_Mb_OrbBased_Tdks

  complex(R64), intent(in), contiguous :: orbs(:, :)
  complex(R64), intent(in), contiguous :: rdm1(:, :)
  real(R64), intent(in) :: time

  complex(R64), allocatable :: rdmSpatial(:, :)
  real(R64) :: seconds, exchangeFraction
  integer(I64) :: count0, count1, countRate
  integer(I32) :: nOS, nPairSolves

  if (.not. referenceSetupQ) &
    error stop "fieldFree.cayley.meanFieldQ: call FieldFree_Setup before building the reference"
  nOS = Orbs_nOrbsInState
  if (size(orbs, 1) .ne. Grid_nPoints .or. size(orbs, 2) .ne. nOS) &
    error stop "FieldFree_BuildMeanField: orbs must be the stored orbitals (nPoints x nOrbsInState)"
  if (size(rdm1, 1) .ne. Method_Mb_OrbBased_nOrbsSum .or. size(rdm1, 2) .ne. Method_Mb_OrbBased_nOrbsSum) &
    error stop "FieldFree_BuildMeanField: rdm1 must be nOrbsSum x nOrbsSum"

  call system_clock(count0, countRate)
  call ReleaseMeanField
  call Method_Mb_OrbBased_FillSpatialRdm1(rdmSpatial, rdm1)

  if (Method_kohnShamQ) then
    ! Kohn-Sham: the local v_KS,s of the state on the xc points; exchange
    ! only with the exact-exchange fraction a of a hybrid, the ACE vectors
    ! carry it (a K0 = sum zeta zeta^H)
    call Method_Mb_OrbBased_Tdks_FillKohnShamPotential(meanFieldKohnSham, meanFieldSpinOfBodyType, orbs, rdm1, time)
    exchangeFraction = SysXc_exactExchangeFraction
    if (exchangeFraction < 0.0_R64) error stop "fieldFree.cayley.meanFieldQ: negative exact-exchange fraction"
    meanFieldExchangeQ = exchangeFraction > 0.0_R64
    nPairSolves = 0
    if (meanFieldExchangeQ) then
      call Method_Mb_OrbBased_FillMeanField(ref, orbs, rdm1, time, nPairSolves_=nPairSolves, exchangeOnlyQ_=.true.)
      ref % zeta(:, :) = sqrt(exchangeFraction) * ref % zeta(:, :)
    end if
  else
    call Method_Mb_OrbBased_FillMeanField(ref, orbs, rdm1, time, nPairSolves_=nPairSolves)
  end if

  FieldFree_meanFieldBuiltQ = .true.

  call system_clock(count1)
  seconds = real(count1 - count0, kind=R64) / real(countRate, kind=R64)
  write (*, '(3A,F10.4,A,F14.10,A,I0,A,F8.2,A)') " exact field-free part: ", ReferenceName(), " built at t =", time, &
    " from the 1-RDM of the state (N =", real(SumDiagonal(rdmSpatial), kind=R64), "): 1 + ", nPairSolves, &
    " potential solves,", seconds, " s"

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module subroutine ApplyMeanField(dOrbs, orbs, h1_)
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

  if (.not. FieldFree_meanFieldBuiltQ) &
    error stop "FieldFree_ApplyMeanField: the frozen mean-field reference is not built &
               &(FieldFree_BuildMeanField)"

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
module subroutine ApplyMeanFieldOnOrb(dOrb, scratch, orb, bt)
  !---------------------------------------------------------------------------
  ! dOrb += U orb of body type bt (AddMeanField) once the reference is built
  ! and FieldFree_meanFieldQ is set (the MCTDHX nonlinear part subtracts
  ! U under the same flag); a no-op otherwise.
  !---------------------------------------------------------------------------
  complex(R64), intent(inout), contiguous :: dOrb(:)
  complex(R64), intent(inout), contiguous :: scratch(:)
  complex(R64), intent(in), contiguous :: orb(:)
  integer(I32), intent(in) :: bt

  if (FieldFree_meanFieldBuiltQ .and. FieldFree_meanFieldQ) call AddMeanField(dOrb, scratch, orb, bt)

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> w = w + U z with the reference of body type bt (Hartree-Fock:
!> Method_Mb_OrbBased_AddMeanField; scratch: work vector). Reentrant; its
!> passes are parallel when called outside a parallel region.
subroutine AddMeanField(w, scratch, z, bt)
  use M_Grid
  use M_SysXc

  complex(R64), intent(inout), contiguous :: w(:)
  complex(R64), intent(inout), contiguous :: scratch(:)
  complex(R64), intent(in), contiguous :: z(:)
  integer(I32), intent(in) :: bt

  complex(R64) :: overlap
  integer(I32) :: k

  if (.not. allocated(meanFieldKohnSham)) then
    call Method_Mb_OrbBased_AddMeanField(w, scratch, ref, z, bt)
    return
  end if

  ! Kohn-Sham: w += v_KS,s z - sum_k zeta_k <zeta_k|z> (zeta scaled by sqrt(a))
  call SysXc_MultiplyWithPotential(scratch, meanFieldKohnSham(:, meanFieldSpinOfBodyType(bt)), z)
  call Axpy(w, (1.0_R64, 0.0_R64), scratch)
  if (.not. meanFieldExchangeQ) return
  do k = ref % zetaStart(bt), ref % zetaEnd(bt)
    overlap = Grid_InnerProduct(ref % zeta(:, k), z)
    call Axpy(w, -overlap, ref % zeta(:, k))
  end do

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Releases the frozen mean-field reference.
subroutine ReleaseMeanField

  FieldFree_meanFieldBuiltQ = .false.
  meanFieldExchangeQ = .true.
  if (allocated(meanFieldKohnSham)) deallocate (meanFieldKohnSham)
  if (allocated(meanFieldSpinOfBodyType)) deallocate (meanFieldSpinOfBodyType)
  ref = T_Method_Mb_OrbBased_MeanField()

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Name of the reference in the log lines.
function ReferenceName() result(res)
  use M_Method, only: Method_kohnShamQ

  character(len=:), allocatable :: res

  if (Method_kohnShamQ) then
    res = "frozen Kohn-Sham reference U = v_KS - a K0"
  else
    res = "frozen mean-field reference U = J0 - K0"
  end if

end function

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
