! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the SingletTr0/Reduce purification backend.
submodule(M_Purification_SingletTr0_Reduce) S_Purification_SingletTr0_Reduce

  implicit none

  !---------------------------------------------------------------------------
  ! Backend configuration (module-scope, set during fabrication)
  !---------------------------------------------------------------------------
  real(R64)    :: threshold       !< Eigenvalue cutoff for defect detection
  integer(I32) :: nSteps          !< Maximum iterations before forced exit
  logical      :: conserveEnergyQ !< If .TRUE., project corrections ⊥ to h²

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Configure the SingletTr0/Reduce backend and bind procedure pointers.
  module subroutine Purification_SingletTr0_Reduce_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Purification

    call Say_Fabricate("purification.singletTr0.reduce")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    nSteps = Json_Get("purification.singletTr0.reduce.nSteps", 40)
    threshold = Json_Get("purification.singletTr0.reduce.threshold", -1e-5_R64)
    conserveEnergyQ = Json_Get("purification.singletTr0.reduce.conserveEnergyQ", .FALSE.)

    Purification_PurifyTwoRdmState => PurifyTwoRdmState

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Purify a SingletTr0-symmetry 2-RDM state by iterative defect reduction.
  !>
  !> @details
  !> Similar to the Singlet/Reduce algorithm, but operates on the trace-free
  !> representation and uses `RdmEnforceTraceM2ab` for the final constraint.
  !> All corrections are accumulated before being applied to rdm2ab.
  !>
  !> @param[inout] twoRdmState  Packed 2-RDM (SingletTr0 layout).
  !> @param[in]    h2           Two-body Hamiltonian for optional energy conservation.
  subroutine PurifyTwoRdmState(twoRdmState, h2)
    use M_TwoRdm_Fermionic_SingletTr0
    use M_Utils_RdmHole
    use M_Utils_RdmAlgebra
    use M_Utils_RdmDefective
    use M_Utils_RdmEnforceTraceM2sym
    use M_Utils_RdmEnforceTraceM2anti
    use M_Utils_RdmProjection
    use M_Utils_RdmEnforceTraceM2ab

    complex(R64), intent(inout), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous     :: h2(:, :, :, :)

    ! Local 2-RDM arrays
    complex(R64), allocatable :: rdm2ab(:, :, :, :)        !< Full α-β 2-RDM
    complex(R64), allocatable :: rdm2singlet(:, :, :, :)   !< Singlet (S=0) block
    complex(R64), allocatable :: rdm2triplet(:, :, :, :)   !< Triplet (S=1) block

    ! Local 2-hole-RDM (Q-RDM) arrays
    complex(R64), allocatable :: qrdm2ab(:, :, :, :)       !< Full α-β Q-RDM
    complex(R64), allocatable :: qrdm2singlet(:, :, :, :)  !< Singlet Q-RDM block
    complex(R64), allocatable :: qrdm2triplet(:, :, :, :)  !< Triplet Q-RDM block

    ! Defective (negative eigenspace) parts, made trace-free
    complex(R64), allocatable :: rdm2singletDefectiveTr0(:, :, :, :)
    complex(R64), allocatable :: rdm2tripletDefectiveTr0(:, :, :, :)
    complex(R64), allocatable :: qrdm2singletDefectiveTr0(:, :, :, :)
    complex(R64), allocatable :: qrdm2tripletDefectiveTr0(:, :, :, :)

    ! Energy conservation workspace
    complex(R64), allocatable :: h2Tr0(:, :, :, :)         !< Trace-free h² (α-β block)
    complex(R64), allocatable :: correction(:, :, :, :)    !< Accumulated correction

    ! Minimum eigenvalues across blocks
    real(R64)    :: minSingletRdm2, minSingletQrdm2
    real(R64)    :: minTripletRdm2, minTripletQrdm2
    real(R64)    :: minEval         !< Overall minimum eigenvalue

    integer(I32) :: i               !< Iteration counter
    integer(I32) :: nNegative       !< Count of negative eigenvalues (diagnostic)
    integer(I32) :: nOa, nOb        !< Orbital dimensions

    !---------------------------------------------------------------------------
    ! Step 1: Unpack the 2-RDM from the compressed state vector
    !---------------------------------------------------------------------------
    call TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)

    !---------------------------------------------------------------------------
    ! Optional: Prepare trace-free h² for energy-conserving projection
    ! Note: For SingletTr0, we extract just the α-β block from h²
    !---------------------------------------------------------------------------
    if (conserveEnergyQ) then
      nOa = size(rdm2ab, 1)
      nOb = nOa
      allocate (h2Tr0(nOa, nOb, nOa, nOb))
      h2Tr0 = h2(1:nOa, nOa + 1:nOa + nOb, 1:nOa, nOa + 1:nOa + nOb)
      call RdmEnforceTraceM2ab_Apply(h2Tr0)
    end if

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
      call RdmAlgebra_FillSingletPart(rdm2singlet, rdm2ab)
      call RdmDefective_FillDefectiveRdm2sym(rdm2singletDefectiveTr0, minSingletRdm2, nNegative, rdm2singlet, threshold)
      call RdmEnforceTraceM2sym_Apply(rdm2singletDefectiveTr0)

      call RdmAlgebra_FillTripletPart(rdm2triplet, rdm2ab)
      call RdmDefective_FillDefectiveRdm2antisym(rdm2tripletDefectiveTr0, minTripletRdm2, nNegative, rdm2triplet, threshold)
      call RdmEnforceTraceM2anti_Apply(rdm2tripletDefectiveTr0)

      correction = rdm2singletDefectiveTr0 + rdm2tripletDefectiveTr0

      !------------------------
      ! 2b. Purify Q² (hole-RDM)
      !------------------------
      call RdmHole_FillQrdm2ab(qrdm2ab, rdm2ab)

      call RdmAlgebra_FillSingletPart(qrdm2singlet, qrdm2ab)
      call RdmDefective_FillDefectiveRdm2sym(qrdm2singletDefectiveTr0, minSingletQrdm2, nNegative, qrdm2singlet, threshold)
      call RdmEnforceTraceM2sym_Apply(qrdm2singletDefectiveTr0)

      call RdmAlgebra_FillTripletPart(qrdm2triplet, qrdm2ab)
      call RdmDefective_FillDefectiveRdm2antisym(qrdm2tripletDefectiveTr0, minTripletQrdm2, nNegative, qrdm2triplet, threshold)
      call RdmEnforceTraceM2anti_Apply(qrdm2tripletDefectiveTr0)

      ! Accumulate all corrections (both D² and Q² defects) before applying
      correction = correction + qrdm2singletDefectiveTr0 + qrdm2tripletDefectiveTr0

      if (conserveEnergyQ) call RdmProjection_Apply(correction, h2Tr0)

      rdm2ab = rdm2ab - correction

      !------------------------
      ! 2c. Convergence check
      !------------------------
      minEval = min(minSingletRdm2, minTripletRdm2, minSingletQrdm2, minTripletQrdm2)

      if (i == 1 .and. minEval < threshold) then
        write (*, *)
        write (*, '(A,E16.6)') " PURIFICATION start with min eval", minEval
      end if

    end do

    write (*, '(A,E16.6, A, I6)') " PURIFICATION ended with min eval", minEval, " at step ", i
    write (*, *)

    !---------------------------------------------------------------------------
    ! Step 3: Repack the purified 2-RDM into the state vector
    !---------------------------------------------------------------------------
    call TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdm2ab(twoRdmState, rdm2ab)

  end subroutine

end submodule
