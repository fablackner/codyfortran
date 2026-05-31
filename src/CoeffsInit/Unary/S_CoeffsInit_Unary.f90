! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Submodule S_CoeffsInit_Unary implements the "Unary" CI coefficient
!> initializer for single-configuration (product-state) initial conditions.
!>
!> Physical Interpretation
!> -----------------------
!> The unary initializer sets the CI expansion to the reference configuration
!> |Φ₀⟩, i.e., all particles occupy their lowest-indexed orbitals according
!> to the Fock-state enumeration in ConfigList. This corresponds to:
!>
!>   |Ψ⟩ = |Φ₀⟩  ⟹  c₁ = 1, c_i = 0 for i > 1
!>
!> This is the natural starting point for:
!>   - Imaginary-time propagation (ground-state search)
!>   - Mean-field TDHF/TDHX calculations
!>   - Initial guess for SCF iterations
!>
!> JSON Configuration
!> ------------------
!>   "coeffsInit": { "unary": { } }
!>
!> No additional parameters are required. The empty block signals intent.
submodule(M_CoeffsInit_Unary) S_CoeffsInit_Unary

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Registers the unary initializer by binding `CoeffsInit_Initialize`
  !> in M_CoeffsInit to the local `Initialize` procedure.
  !>
  !> @note No setup phase is needed for this initializer; the default
  !>       no-op `CoeffsInit_Setup` remains in place.
  module subroutine CoeffsInit_Unary_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_CoeffsInit

    call Say_Fabricate("coeffsInit.unary")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    CoeffsInit_Initialize => Initialize

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Populates the CI coefficient vector with the reference configuration.
  !>
  !> Sets c(1) = 1 (normalized) and c(i) = 0 for all i > 1, corresponding
  !> to a pure product state |Φ₀⟩ in the Fock basis.
  !>
  !> @param[out] coeffs  Complex CI coefficient vector; dimension must match
  !>                     the CI space size from ConfigList.
  subroutine Initialize(coeffs)
    use M_Utils_DataStorage

    complex(R64), intent(out), contiguous  :: coeffs(:)

    coeffs(:) = 0.0_R64
    coeffs(1) = 1.0_R64

  end subroutine

end submodule
