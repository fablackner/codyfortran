! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the Singlet/Reduce purification backend.
submodule(M_Purification_Singlet_Reduce) S_Purification_Singlet_Reduce

  implicit none

  !---------------------------------------------------------------------------
  ! Backend configuration (module-scope, set during fabrication)
  !---------------------------------------------------------------------------
  real(R64)    :: threshold       !< Eigenvalue cutoff for defect detection
  integer(I32) :: nSteps          !< Maximum iterations before forced exit
  logical      :: conserveEnergyQ !< If .TRUE., project corrections ⊥ to h²

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Configure the Singlet/Reduce backend and bind procedure pointers.
  module subroutine Purification_Singlet_Reduce_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Purification

    call Say_Fabricate("purification.singlet.reduce")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    nSteps = Json_Get("purification.singlet.reduce.nSteps", 40)
    threshold = Json_Get("purification.singlet.reduce.threshold", -1e-5_R64)
    conserveEnergyQ = Json_Get("purification.singlet.reduce.conserveEnergyQ", .FALSE.)
    Purification_usesH2Q = conserveEnergyQ

    Purification_PurifyTwoRdmState => PurifyTwoRdmState

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Purify a singlet-symmetry 2-RDM state by iterative defect reduction.
  !>
  !> @details
  !> This is the Mazziotti/Alcoba-style iterative purification for singlet
  !> systems: unpack the (up,down) block, reduce it
  !> (Purification_Singlet_Reduce_PurifyRdm2ab), pack it.
  !>
  !> @param[inout] twoRdmState  Packed 2-RDM (Singlet layout: α-β block).
  !> @param[in]    h2           Two-body Hamiltonian for optional energy conservation
  !>                            (full spin-orbital tensor, see M_Purification).
  subroutine PurifyTwoRdmState(twoRdmState, h2)
    use M_TwoRdm_Fermionic_Singlet

    complex(R64), intent(inout), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous     :: h2(:, :, :, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)        !< Full α-β 2-RDM
    complex(R64), allocatable :: h2Tr0(:, :, :, :)         !< Trace-free energy direction

    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)

    if (conserveEnergyQ) then
      call Purification_Singlet_Reduce_FillEnergyDirection(h2Tr0, h2, size(rdm2ab, 1))
      call Purification_Singlet_Reduce_PurifyRdm2ab(rdm2ab, nSteps, threshold, h2Tr0_=h2Tr0)
    else
      call Purification_Singlet_Reduce_PurifyRdm2ab(rdm2ab, nSteps, threshold)
    end if

    call TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab(twoRdmState, rdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Purification_Singlet_Reduce_PurifyRdm2ab(rdm2ab, nSteps, threshold, h2Tr0_)
    use M_Method_Mb
    use M_Utils_RdmHole
    use M_Utils_RdmAlgebra
    use M_Utils_RdmDefective
    use M_Utils_RdmEnforceTraceM2sym
    use M_Utils_RdmEnforceTraceM2anti
    use M_Utils_RdmProjection

    complex(R64), intent(inout), contiguous :: rdm2ab(:, :, :, :)
    integer(I32), intent(in) :: nSteps
    real(R64), intent(in) :: threshold
    complex(R64), intent(in), contiguous, optional :: h2Tr0_(:, :, :, :)

    ! Spin-adapted blocks of D² and of the 2-hole-RDM Q²
    complex(R64), allocatable :: rdm2singlet(:, :, :, :), rdm2triplet(:, :, :, :)
    complex(R64), allocatable :: qrdm2ab(:, :, :, :), qrdm2singlet(:, :, :, :), qrdm2triplet(:, :, :, :)

    ! Defective (negative eigenspace) parts, made trace-free
    complex(R64), allocatable :: rdm2singletDefectiveTr0(:, :, :, :)
    complex(R64), allocatable :: rdm2tripletDefectiveTr0(:, :, :, :)
    complex(R64), allocatable :: qrdm2singletDefectiveTr0(:, :, :, :)
    complex(R64), allocatable :: qrdm2tripletDefectiveTr0(:, :, :, :)
    complex(R64), allocatable :: correction(:, :, :, :)    !< Correction of one step

    ! Minimum eigenvalues across blocks
    real(R64)    :: minSingletRdm2, minSingletQrdm2
    real(R64)    :: minTripletRdm2, minTripletQrdm2
    real(R64)    :: minEval         !< Overall minimum eigenvalue

    integer(I32) :: i               !< Iteration counter
    integer(I32) :: nNegative       !< Count of negative eigenvalues (diagnostic)

    !---------------------------------------------------------------------------
    ! Step 2: Iterative purification loop
    !---------------------------------------------------------------------------
    minEval = -1.0_R64
    i = 0

    do while (minEval < threshold .and. i < nSteps)
      i = i + 1

      !------------------------
      ! 2a. Purify D² (2-RDM)
      !------------------------
      ! Decompose into spin-adapted blocks
      call RdmAlgebra_FillSingletPart(rdm2singlet, rdm2ab)
      call RdmAlgebra_FillTripletPart(rdm2triplet, rdm2ab)

      ! Extract defective (negative eigenvalue) parts and enforce trace-free
      call RdmDefective_FillDefectiveRdm2sym(rdm2singletDefectiveTr0, minSingletRdm2, nNegative, rdm2singlet, threshold)
      call RdmEnforceTraceM2sym_Apply(rdm2singletDefectiveTr0)

      call RdmDefective_FillDefectiveRdm2antisym(rdm2tripletDefectiveTr0, minTripletRdm2, nNegative, rdm2triplet, threshold)
      call RdmEnforceTraceM2anti_Apply(rdm2tripletDefectiveTr0)

      ! Accumulate and apply correction to D²
      correction = rdm2singletDefectiveTr0 + rdm2tripletDefectiveTr0
      if (present(h2Tr0_)) call RdmProjection_Apply(correction, h2Tr0_)
      rdm2ab = rdm2ab - correction

      !------------------------
      ! 2b. Purify Q² (hole-RDM)
      !------------------------
      ! Compute Q-RDM from updated D²
      call RdmHole_FillQrdm2ab(qrdm2ab, rdm2ab, Method_Mb_nBodies(1), Method_Mb_nBodies(2))

      ! Decompose into spin-adapted blocks
      call RdmAlgebra_FillSingletPart(qrdm2singlet, qrdm2ab)
      call RdmAlgebra_FillTripletPart(qrdm2triplet, qrdm2ab)

      ! Extract defective parts and enforce trace-free
      call RdmDefective_FillDefectiveRdm2sym(qrdm2singletDefectiveTr0, minSingletQrdm2, nNegative, qrdm2singlet, threshold)
      call RdmEnforceTraceM2sym_Apply(qrdm2singletDefectiveTr0)

      call RdmDefective_FillDefectiveRdm2antisym(qrdm2tripletDefectiveTr0, minTripletQrdm2, nNegative, qrdm2triplet, threshold)
      call RdmEnforceTraceM2anti_Apply(qrdm2tripletDefectiveTr0)

      ! Note: Q-RDM correction applied to D² (not Q²) due to the mapping
      correction = qrdm2singletDefectiveTr0 + qrdm2tripletDefectiveTr0
      if (present(h2Tr0_)) call RdmProjection_Apply(correction, h2Tr0_)
      rdm2ab = rdm2ab - correction

      !------------------------
      ! 2c. Convergence check
      !------------------------
      minEval = min(minSingletRdm2, minTripletRdm2, minSingletQrdm2, minTripletQrdm2)

      if (i .eq. 1 .and. minEval < threshold) then
        write (*, *)
        write (*, '(A,E16.6)') " PURIFICATION start with min eval", minEval
      end if

    end do

    write (*, '(A,E16.6, A, I6)') " PURIFICATION ended with min eval", minEval, " at step ", i
    write (*, *)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Purification_Singlet_Reduce_FillEnergyDirection(h2Tr0, h2, nOa)
    use M_Utils_RdmAlgebra
    use M_Utils_RdmEnforceTraceM2sym
    use M_Utils_RdmEnforceTraceM2anti

    complex(R64), intent(out), allocatable :: h2Tr0(:, :, :, :)
    complex(R64), intent(in), contiguous :: h2(:, :, :, :)
    integer(I32), intent(in) :: nOa

    complex(R64), allocatable :: h2ab(:, :, :, :), tmp(:, :, :, :)

    ! The (up,down) block of the spin-orbital h2
    if (size(h2, 1) .ne. 2 * nOa) error stop "Purification_Singlet_Reduce: h2 must be the spin-orbital tensor (2 nOa)^4"
    h2ab = h2(1:nOa, nOa + 1:2 * nOa, 1:nOa, nOa + 1:2 * nOa)

    call RdmAlgebra_FillSingletPart(tmp, h2ab)
    call RdmEnforceTraceM2sym_Apply(tmp)
    h2Tr0 = tmp
    deallocate (tmp)

    call RdmAlgebra_FillTripletPart(tmp, h2ab)
    call RdmEnforceTraceM2anti_Apply(tmp)
    h2Tr0 = h2Tr0 + 3.0_R64 * tmp

  end subroutine

end submodule
