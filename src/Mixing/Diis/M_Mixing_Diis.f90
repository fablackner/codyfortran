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
!>   1. Stores the residual r = xRaw − xMixed and the damped candidate
!>      y = xMixed + α·r in a circular buffer of capacity `nHistory`, and the
!>      new row Re⟨r_i|r⟩ of the Pulay matrix B (kept across calls: one
!>      matrix-vector product per step)
!>   2. Solves the Pulay problem: minimize ‖Σ_i c_i·r_i‖² subject to
!>      Σ_i c_i = 1, via the bordered system
!>
!>        [ B   1 ] [c]   [0]
!>        [ 1ᵀ  0 ] [λ] = [1]
!>
!>      with the pseudo-inverse of its symmetric eigen-decomposition
!>      (eigenvalues below 1e-12 relative dropped), which stays
!>      well-behaved for (nearly) linearly dependent residual histories
!>   3. Takes the damped step along the extrapolated residual:
!>
!>        xMixed ← Σ_i c_i·(xMixed_i + α·r_i) = Σ_i c_i·y_i
!>
!> With a single stored pair this is exactly linear mixing
!> xMixed ← (1−α)·xMixed + α·xRaw.
!>
!> Safeguards (all keep only the newest pair, i.e. take a plain damped step):
!>   - **start gate**: plain steps while the residual norm is above
!>     `startThreshold` × the first residual norm (the linear model holds
!>     only near convergence); `startThreshold` 0 switches the gate off
!>     (extrapolation from the second step on)
!>   - **growth**: the residual norm grows by more than 3× in one step
!>   - **coefficients**: Σ|c_i| exceeds 1e4 (near-dependent residuals)
!>
!> **JSON Configuration:**
!> @code{.json}
!> "mixing": { "diis": { "nHistory": 8, "startThreshold": 0.1, "alpha": 1.0 } }
!> @endcode
!>
!> | Parameter        | Type | Default | Description                               |
!> |------------------|------|---------|-------------------------------------------|
!> | `nHistory`       | int  | 8       | Maximum stored pairs                      |
!> | `startThreshold` | real | 0.1     | Relative residual drop below which the    |
!> |                  |      |         | extrapolation acts (0: no gate)           |
!> | `alpha`          | real | 1.0     | Damping parameter (0 < α ≤ 1)            |
!> | `printLevel`     | int  | 0       | 1: one line per step (residual norm,      |
!> |                  |      |         | stored pairs, gate and restarts)          |
!>
!> Memory: 2·nHistory vectors of the mixed length.
!>
!> @note DIIS assumes that the mixed quantity determines the next iterate:
!>       the fixed-point map must be a function of it alone. Orbitals and
!>       the Kohn-Sham potential of a local functional are such iterates; the
!>       Hartree potential of Hartree-Fock (or of a hybrid, or the MCSCF
!>       mean field) is not, because the exchange is rebuilt from the
!>       orbitals and lags the mixed potential (A_Mixing.md, "What to mix").
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
