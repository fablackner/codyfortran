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
  !> This is the core algorithm implementing the Mazziotti/Alcoba-style iterative
  !> purification for singlet systems. The 2-RDM is decomposed into spin-adapted
  !> blocks and both D²- and Q²-conditions are enforced simultaneously.
  !>
  !> @param[inout] twoRdmState  Packed 2-RDM (Singlet layout: α-β block).
  !> @param[in]    h2           Two-body Hamiltonian for optional energy conservation
  !>                            (full spin-orbital tensor, see M_Purification;
  !>                            the (up,down) block alone is accepted as well).
  subroutine PurifyTwoRdmState(twoRdmState, h2)
    use M_TwoRdm_Fermionic_Singlet
    use M_Utils_RdmHole
    use M_Utils_RdmAlgebra
    use M_Utils_RdmDefective
    use M_Utils_RdmEnforceTraceM2sym
    use M_Utils_RdmEnforceTraceM2anti
    use M_Utils_RdmProjection

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
    complex(R64), allocatable :: h2ab(:, :, :, :)          !< (up,down) block of h²
    complex(R64), allocatable :: h2Tr0(:, :, :, :)         !< Trace-free energy functional
    complex(R64), allocatable :: correction(:, :, :, :)    !< Accumulated correction
    complex(R64), allocatable :: tmp(:, :, :, :)           !< Temporary workspace

    ! Minimum eigenvalues across blocks
    real(R64)    :: minSingletRdm2, minSingletQrdm2
    real(R64)    :: minTripletRdm2, minTripletQrdm2
    real(R64)    :: minEval         !< Overall minimum eigenvalue

    integer(I32) :: i               !< Iteration counter
    integer(I32) :: nNegative       !< Count of negative eigenvalues (diagnostic)
    integer(I32) :: nOa             !< Number of spin-up orbitals

    !---------------------------------------------------------------------------
    ! Step 1: Unpack the 2-RDM from the compressed state vector
    !---------------------------------------------------------------------------
    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)

    !---------------------------------------------------------------------------
    ! Optional: Prepare trace-free h² for energy-conserving projection
    !---------------------------------------------------------------------------
    if (conserveEnergyQ) then
      ! The (up,down) block of h2; a full spin-orbital h2 is reduced to it
      nOa = size(rdm2ab, 1)
      if (size(h2, 1) .eq. nOa) then
        allocate (h2ab, source=h2)
      else if (size(h2, 1) .eq. 2 * nOa) then
        h2ab = h2(1:nOa, nOa + 1:2 * nOa, 1:nOa, nOa + 1:2 * nOa)
      else
        error stop "Purification_Singlet_Reduce: h2 matches neither the (up,down) block nor the spin-orbital basis"
      end if

      ! For a singlet the up-up block follows from the (up,down) block,
      ! D_aa = D_ab - P12 D_ab, so the two-body energy as a functional of
      ! D_ab is <4 h2 - 2 P12 h2 | D_ab> = <2 h2_singlet + 6 h2_triplet | D_ab>.
      ! Corrections orthogonal to (trace-free) h2_singlet + 3 h2_triplet
      ! therefore conserve the energy (the 1-RDM is unchanged anyway).
      allocate (h2Tr0, mold=h2ab)
      h2Tr0 = 0.0_R64

      call RdmAlgebra_FillSingletPart(tmp, h2ab)
      call RdmEnforceTraceM2sym_Apply(tmp)
      h2Tr0 = h2Tr0 + tmp
      deallocate (tmp)

      call RdmAlgebra_FillTripletPart(tmp, h2ab)
      call RdmEnforceTraceM2anti_Apply(tmp)
      h2Tr0 = h2Tr0 + 3.0_R64 * tmp

      if (allocated(tmp)) deallocate (tmp)
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
      if (conserveEnergyQ) call RdmProjection_Apply(correction, h2Tr0)
      rdm2ab = rdm2ab - correction

      !------------------------
      ! 2b. Purify Q² (hole-RDM)
      !------------------------
      ! Compute Q-RDM from updated D²
      call RdmHole_FillQrdm2ab(qrdm2ab, rdm2ab)

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
    call TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab(twoRdmState, rdm2ab)

  end subroutine

end submodule
