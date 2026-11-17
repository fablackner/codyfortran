! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Submodule S_CoeffsInit_Excited implements the "Excited" CI coefficient
!> initializer that constructs states via creation/annihilation operators.
!>
!> Physical Interpretation
!> -----------------------
!> This initializer builds an excited configuration by applying a sequence
!> of second-quantized operators to the reference state |Φ₀⟩:
!>
!>   |Ψ⟩ = â†_{p₁} â†_{p₂} ... â_{q₁} â_{q₂} ... |Φ₀⟩
!>
!> where `creates` specifies target orbitals p_i and `destroys` specifies
!> source orbitals q_i. The operation is applied independently to each
!> body type (e.g., spin-up and spin-down for fermions).
!>
!> Use Cases
!> ---------
!>   - Initialize specific excited configurations for dynamics
!>   - Prepare particle-hole excitations for response calculations
!>   - Build spin-flip or charge-transfer states
!>
!> JSON Configuration
!> ------------------
!>   "coeffsInit": {
!>     "excited": {
!>       "creates":   [3, 4],    // orbital indices to create particles in
!>       "destroys":  [1, 2],    // orbital indices to annihilate from
!>       "bodyType1": 1,         // first body type (default: 1)
!>       "bodyType2": 2          // second body type (default: 2)
!>     }
!>   }
!>
!> The excitation â†_{creates} â_{destroys} is applied to bodyType1 and,
!> if different, to bodyType2 (default: body type 2, or 1 for a single body
!> type), and the result is normalized.
submodule(M_CoeffsInit_Excited) S_CoeffsInit_Excited

  implicit none

  !> Target orbital indices for particle creation operators â†_p.
  !> Read from JSON: "coeffsInit.excited.creates"
  integer(I32), allocatable :: creates(:)

  !> Source orbital indices for annihilation operators â_q.
  !> Read from JSON: "coeffsInit.excited.destroys"
  integer(I32), allocatable :: destroys(:)

  !> First body type index (e.g., spin-up electrons, species A).
  !> Default: 1
  integer(I32) :: bodyType1 = 1

  !> Second body type index (e.g., spin-down electrons, species B).
  !> Default: 2 (1 for a single body type)
  integer(I32) :: bodyType2 = 2

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Registers the excited-state initializer and reads JSON parameters.
  !>
  !> Parses:
  !>   - coeffsInit.excited.creates   → creation operator targets
  !>   - coeffsInit.excited.destroys  → annihilation operator sources
  !>   - coeffsInit.excited.bodyType1 → first body type (default 1)
  !>   - coeffsInit.excited.bodyType2 → second body type (default 2, or 1
  !>     for a single body type)
  module subroutine CoeffsInit_Excited_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Method_Mb
    use M_CoeffsInit

    call Say_Fabricate("coeffsInit.excited")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    creates = Json_Get("coeffsInit.excited.creates", [1])
    destroys = Json_Get("coeffsInit.excited.destroys", [1])
    bodyType1 = Json_Get("coeffsInit.excited.bodyType1", 1)
    bodyType2 = Json_Get("coeffsInit.excited.bodyType2", min(2, Method_Mb_nBodyTypes))
    if (bodyType1 < 1 .or. bodyType1 > Method_Mb_nBodyTypes .or. bodyType2 < 1 .or. bodyType2 > Method_Mb_nBodyTypes) &
      error stop "coeffsInit.excited: bodyType1 and bodyType2 must be body types of method.mb"

    CoeffsInit_Initialize => Initialize

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Constructs the excited-state CI vector by applying excitation operators.
  !>
  !> Algorithm:
  !>   1. Initialize to reference state: c(1) = 1, c(i>1) = 0
  !>   2. Apply â†_{creates} â_{destroys} for bodyType1
  !>   3. Apply â†_{creates} â_{destroys} for bodyType2 if it differs
  !>   4. Normalize (bosonic ladder factors); stop if the excitation
  !>      annihilated the reference configuration
  !>
  !> The `Coeffs_ApplyExcitation` routine handles the Fock-space algebra,
  !> including fermionic sign factors and configuration lookup.
  !>
  !> @param[out] coeffs  Complex CI coefficient vector to be filled.
  subroutine Initialize(coeffs)
    use M_Coeffs

    complex(R64), intent(out), contiguous :: coeffs(:)

    ! Start from the reference configuration |Φ₀⟩
    coeffs(:) = 0.0_R64
    coeffs(1) = 1.0_R64

    ! Apply excitation operators for both body types (once for one)
    call Coeffs_ApplyExcitation(coeffs, creates, destroys, bodyType1)
    if (bodyType2 .ne. bodyType1) call Coeffs_ApplyExcitation(coeffs, creates, destroys, bodyType2)

    if (.not. (sum(abs(coeffs)**2) > 0.0_R64)) &
      error stop "coeffsInit.excited: the excitation annihilates the reference configuration (zero CI vector)"
    call Coeffs_Normalize(coeffs)

  end subroutine

end submodule
