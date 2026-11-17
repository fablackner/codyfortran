! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Cayley variant of the exact field-free part (grids whose kinetic backend
!> binds SysKinetic_SolveShifted, e.g. grid.fem): every orbital is propagated
!> by the Crank-Nicolson step (1 + i dt h0/2)^-1 (1 - i dt h0/2), solved by
!> right-preconditioned GMRES in the grid metric (M_Utils_Gmres). Unitary and
!> second order; real time only; no phi functions and no spectral cut.
!>
!> JSON ("fieldFree.cayley"):
!>   - tolerance:   relative residual of the GMRES solve (default 1e-12)
!>   - krylovDim:   restart length (default 30)
!>   - maxRestarts: restarts (default 20)
!>   - meanFieldQ:  frozen mean-field reference U (default false): the step
!>                  propagates h0 + U, U = J0 - K0 of the state at the first
!>                  propagation (U = v_KS[rho0] - a K0 for a Kohn-Sham
!>                  method), and the method subtracts U in its nonlinear part
!>                  (FieldFree_BuildMeanField, FieldFree_ApplyMeanField;
!>                  S_FieldFree_Cayley_MeanField). Fermions; a
!>                  body-type-independent interaction.
!>
!> Requires a time-independent kinetic operator and a time- and
!> body-type-independent potential (a laser belongs to sysGauge).
module M_FieldFree_Cayley
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Reads "fieldFree.cayley" and binds FieldFree_Setup, FieldFree_Propagate
    !> and the mean-field hooks FieldFree_BuildMeanField and
    !> FieldFree_ApplyMeanField.
    module subroutine FieldFree_Cayley_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Diagnostics of the Cayley step: number of GMRES solves (one per orbital
  !> and exact substep), their total and largest iteration counts
  !> (preconditioned matrix-vector products, all restarts) since FieldFree_Setup.
  integer(I64) :: FieldFree_Cayley_nSolves = 0
  integer(I64) :: FieldFree_Cayley_nIterations = 0
  integer(I32) :: FieldFree_Cayley_maxIterations = 0

end module
