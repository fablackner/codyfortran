! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Generic spin‑orbital fermionic backend for 2‑RDM operations.
!>
!> @details
!> This backend treats all spin‑orbitals explicitly (unrestricted) and applies
!> standard fermionic antisymmetry and contraction rules. No spin symmetry is
!> exploited, making it suitable for broken-symmetry states or benchmarking.
!>
!> ## Storage Layout
!>
!> The linearized 2-RDM uses the full nOrbs⁴ tensor in column-major order:
!> ```
!> twoRdmState(count) = rdm2(i₁, i₂, j₁, j₂)
!> count = i₁ + nO*(i₂-1) + nO²*(j₁-1) + nO³*(j₂-1)
!> ```
!>
!> ## Index Convention
!>
!> The 2-RDM tensor follows the physicist's convention:
!>   D²(i₁, i₂; j₁, j₂) = ⟨Ψ| a†_{j₁} a†_{j₂} a_{i₂} a_{i₁} |Ψ⟩
!>
!> Antisymmetry: D²(i₁,i₂;j₁,j₂) = -D²(i₂,i₁;j₁,j₂) = -D²(i₁,i₂;j₂,j₁)
!>
!> ## Commutator Algebra
!>
!> One-body: [ĥ₁, D²]_{i₁i₂;j₁j₂} = Σ_k (h_{i₁k}D²_{ki₂;j₁j₂} + h_{i₂k}D²_{i₁k;j₁j₂}
!>                                       - h_{kj₁}D²_{i₁i₂;kj₂} - h_{kj₂}D²_{i₁i₂;j₁k})
!>
!> Two-body: [Ŵ₁₂, D²] = 2·Σ_{k₁k₂} (W_{i₁i₂;k₁k₂}D²_{k₁k₂;j₁j₂} - W_{k₁k₂;j₁j₂}D²_{i₁i₂;k₁k₂})
module M_TwoRdm_Fermionic_Generic
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Activate the generic fermionic backend and wire procedure pointers.
    !>
    !> Responsibilities:
    !> - Assign `M_TwoRdm` procedure pointers to the generic implementations.
    !> - Define/record the serialization of the linearized 2‑RDM used by this backend.
    !> - Prepare any precomputed index maps or symmetry factors required at runtime.
    module subroutine TwoRdm_Fermionic_Generic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! local data
  !=============================================================================

end module
