! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief DIIS (Pulay/Anderson) mixing backend.
!>
!> @details
!> Accelerates fixed-point iterations by extrapolating over a history of
!> (iterate, residual) pairs. Each `Mixing_Mix` call:
!>
!>   1. Stores the pair (xMixed, r) with residual r = xRaw − xMixed in a circular
!>      buffer of capacity `nHistory`
!>   2. Solves the Pulay problem: minimize ‖Σ_i c_i·r_i‖² subject to
!>      Σ_i c_i = 1, via the bordered linear system
!>
!>        [ B   1 ] [c]   [0]        B_ij = ⟨r_i|r_j⟩
!>        [ 1ᵀ  0 ] [λ] = [1]
!>
!>      solved with an SVD pseudo-inverse (`LapackLib_Svd`), which stays
!>      well-behaved for (nearly) linearly dependent residual histories
!>   3. Takes a damped step along the extrapolated residual:
!>
!>        xMixed ← Σ_i c_i·(xMixed_i + α·r_i)
!>
!> With a single stored pair this reduces exactly to linear mixing
!> xMixed ← (1−α)·xMixed + α·xRaw, so the first iteration behaves like the linear
!> backend and acceleration builds up as history accumulates.
!>
!> Safeguards (standard for production SCF codes): extrapolation is gated
!> until the residual norm has dropped below `startThreshold` × (first
!> residual norm) — before that, and whenever the residual grows sharply or
!> the Pulay coefficients blow up, the history is restarted and a plain
!> damped step is taken.
!>
!> **JSON Configuration:**
!> @code{.json}
!> "mixing": { "diis": { "nHistory": 8, "startThreshold": 0.1, "alpha": 1.0 } }
!> @endcode
!>
!> | Parameter        | Type | Default | Description                               |
!> |------------------|------|---------|-------------------------------------------|
!> | `nHistory`       | int  | 8       | Maximum stored (iterate, residual) pairs  |
!> | `startThreshold` | real | 0.1     | Relative residual drop that enables       |
!> |                  |      |         | extrapolation                             |
!> | `alpha`          | real | 1.0     | Damping parameter (0 < α ≤ 1)            |
!>
!> @note DIIS assumes the mixed quantity lives in a linear space and that the
!>       fixed-point map is smooth in it. Densities and potentials satisfy
!>       this; orbital sets only do so after gauge alignment
!>       (`Orbs_AlignOnReference`) and may still converge more slowly
!>       than potential-target DIIS because of the orthonormalization step.
module M_Mixing_Diis
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface Mixing_Diis_Fabricate
    !> @brief Read JSON parameters and bind the DIIS implementation.
    module subroutine Mixing_Diis_Fabricate()
    end subroutine
  end interface

end module
