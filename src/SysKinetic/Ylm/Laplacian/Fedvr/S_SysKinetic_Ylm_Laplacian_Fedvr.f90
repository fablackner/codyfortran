! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Ylm_Laplacian_Fedvr.f90
!> @brief Implementation of FEDVR radial kinetic operator for Ylm channels.
submodule(M_SysKinetic_Ylm_Laplacian_Fedvr) S_SysKinetic_Ylm_Laplacian_Fedvr
  !-----------------------------------------------------------------------------
  ! Besides the radial kinetic operator, this backend binds the optional exact
  ! shifted solve SysKinetic_SolveShifted, (T + shift)^-1. T is block diagonal
  ! in (l,m) and independent of m; in the metric-symmetrized radial basis
  ! (w^{1/2} f, w the radial metric weights) every block of body type 1 is
  ! real symmetric, ts_l = U_l E_l U_l^T, so that
  !   (T + shift)^-1 f_lm = w^{-1/2} U_l (E_l m_1/m_bt + shift)^-1 U_l^T w^{1/2} f_lm
  ! (T is proportional to 1/m). The blocks are built (probing with radial unit
  ! vectors) and diagonalized (LapackLib_DiagonalizeSymStack, as the exact
  ! field-free part of Method does for h0) at the first call, so that runs
  ! without a consumer pay nothing.
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Eigenvectors U_l of the symmetrized kinetic blocks (body type 1),
  !> (nRadial, nRadial, 0:lmax).
  real(R64), allocatable :: kineticEvecs(:, :, :)

  !> Eigenvalues E_l of the symmetrized kinetic blocks (body type 1),
  !> (nRadial, 0:lmax).
  real(R64), allocatable :: kineticEvals(:, :)

  !> Square roots of the radial metric weights.
  real(R64), allocatable :: sqrtWeights(:)

  !> Whether the kinetic blocks are diagonalized (SetupSolveShifted done).
  logical :: solveShiftedReadyQ = .false.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind FEDVR radial operator, validate grid requirements.
  module subroutine SysKinetic_Ylm_Laplacian_Fedvr_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Ylm

    call Say_Fabricate("sysKinetic.ylm.laplacian.fedvr")

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    SysKinetic_Ylm_MultiplyWithRadialKineticOp => MultiplyWithRadialKineticOp
    SysKinetic_SolveShifted => SolveShifted

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.ylm.fedvr")) then
      error stop "grid.ylm.fedvr is required for sysKinetic.ylm.laplacian.fedvr"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply the radial kinetic operator using FEDVR and g(r) = r·f(r) transformation.
  !>
  !> @details
  !> Algorithm (same as FinDiff, but uses FEDVR derivative context):
  !>   1. Transform: g(r) = r · f(r)
  !>   2. Compute d²g/dr² using FEDVR derivative matrices
  !>   3. Inverse transform: result = (1/r) d²g/dr²
  !>   4. Add centrifugal term: − l(l+1)/r² · f(r)
  !>   5. Scale by −1/(2m)
  subroutine MultiplyWithRadialKineticOp(dOrbLm, orbLm, l, m, time, bt_)
    use M_Utils_Fedvr
    use M_Utils_UnusedVariables
    use M_Utils_DerivativeFedvr
    use M_Grid
    use M_Grid_Ylm
    use M_Grid_Ylm_Fedvr, fedvrCtx => Grid_Ylm_Fedvr_fedvrCtx, derivativeCtx => Grid_Ylm_Fedvr_derivativeCtx
    use M_SysKinetic_Ylm_Laplacian

    complex(R64), intent(out) :: dOrbLm(:)
    complex(R64), intent(in)  :: orbLm(:)
    integer(I32), intent(in)  :: l
    integer(I32), intent(in)  :: m
    real(R64), intent(in)     :: time
    integer(I32), intent(in), optional :: bt_

    complex(R64), allocatable :: gRadial(:)
    integer(I32) :: bt

    if (.false.) call UnusedVariables_Mark(m, time)

    bt = 1
    if (present(bt_)) bt = bt_

    allocate (gRadial(Grid_Ylm_nRadial))
    gRadial = 0.0_R64

    ! Transform: g(r) = r · f(r)
    gRadial(:) = Grid_Ylm_radialPoints(:) * orbLm(:)

    ! Compute d²g/dr² using FEDVR method
    call DerivativeFedvr_Do2ndDerivative(dOrbLm, gRadial, derivativeCtx, fedvrCtx)

    ! Inverse transform: result = (1/r) · d²g/dr²
    dOrbLm(:) = dOrbLm(:) / Grid_Ylm_radialPoints(:)

    ! Add centrifugal term: −l(l+1)/r² · f(r)
    if (l > 0) then
      dOrbLm(:) = dOrbLm(:) - (l * (l + 1)) * orbLm(:) / Grid_Ylm_radialPoints(:)**2
    end if

    ! Scale by −1/(2m)
    dOrbLm(:) = -0.5_R64 * dOrbLm(:) / SysKinetic_Ylm_Laplacian_bodyMass(bt)

    deallocate (gRadial)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Exact shifted solve dOrb = (T̂ + shift)⁻¹·orb (SysKinetic_SolveShifted).
  !>
  !> @details
  !> Channel by channel with the eigen-decomposition of the radial kinetic
  !> block of l (see the submodule header); diagonalizes the blocks at the
  !> first call. T̂ + shift must be invertible: for a real shift,
  !> shift > −(smallest eigenvalue of T̂) (> 0 always works).
  subroutine SolveShifted(dOrb, orb, shift, bt_)
    use M_Grid_Ylm
    use M_SysKinetic_Ylm_Laplacian

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous  :: orb(:)
    complex(R64), intent(in)              :: shift
    integer(I32), intent(in), optional    :: bt_

    complex(R64), allocatable :: fLm(:), denominators(:)
    real(R64) :: massRatio
    integer(I32) :: bt, l, m
    logical :: readyQ

    bt = 1
    if (present(bt_)) bt = bt_

    !$omp atomic read
    readyQ = solveShiftedReadyQ
    if (.not. readyQ) then
      !$omp critical (SysKinetic_Ylm_Laplacian_Fedvr_SetupSolveShifted)
      if (.not. solveShiftedReadyQ) call SetupSolveShifted
      !$omp end critical (SysKinetic_Ylm_Laplacian_Fedvr_SetupSolveShifted)
    end if

    massRatio = SysKinetic_Ylm_Laplacian_bodyMass(1) / SysKinetic_Ylm_Laplacian_bodyMass(bt)

    allocate (fLm(Grid_Ylm_nRadial), denominators(Grid_Ylm_nRadial))
    dOrb(:) = (0.0_R64, 0.0_R64)

    do l = 0, Grid_Ylm_lmax
      denominators(:) = 1.0_R64 / (massRatio * kineticEvals(:, l) + shift)
      do m = -min(l, Grid_Ylm_mmax), min(l, Grid_Ylm_mmax)
        call Grid_Ylm_GetLmComponent(fLm, l, m, orb)
        ! w^{1/2} f -> U^T (as the row vector product f^T U) -> /(E + shift)
        ! -> U -> w^{-1/2}
        fLm(:) = sqrtWeights(:) * fLm(:)
        fLm(:) = denominators(:) * matmul(fLm, kineticEvecs(:, :, l))
        fLm(:) = matmul(kineticEvecs(:, :, l), fLm) / sqrtWeights(:)
        call Grid_Ylm_SetLmComponent(dOrb, l, m, fLm)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Builds and diagonalizes the metric-symmetrized radial kinetic
  !> blocks of body type 1 for l = 0..lmax (once, from SolveShifted).
  subroutine SetupSolveShifted
    use M_Utils_LapackLib
    use M_Grid_Ylm

    real(R64), parameter :: tolerance = 1.0e-10_R64

    complex(R64), allocatable :: unitRad(:), tLm(:)
    real(R64) :: devSym
    integer(I32) :: nRad, lmax, i, l

    if (Grid_Ylm_cProductQ) error stop "SysKinetic_SolveShifted (ylm.laplacian.fedvr) requires a Hermitian (non-ECS) grid"

    nRad = Grid_Ylm_nRadial
    lmax = Grid_Ylm_lmax

    allocate (sqrtWeights(nRad))
    sqrtWeights(:) = sqrt(Grid_Ylm_radialWeights(:))

    allocate (unitRad(nRad), tLm(nRad))
    allocate (kineticEvecs(nRad, nRad, 0:lmax), kineticEvals(nRad, 0:lmax))

    ! Column i of the block of l in metric-symmetrized form
    do l = 0, lmax
      do i = 1, nRad
        unitRad(:) = (0.0_R64, 0.0_R64)
        unitRad(i) = (1.0_R64, 0.0_R64)
        call MultiplyWithRadialKineticOp(tLm, unitRad, l, 0, 0.0_R64, bt_=1)
        kineticEvecs(:, i, l) = sqrtWeights(:) * real(tLm(:), kind=R64) / sqrtWeights(i)
      end do
    end do

    call LapackLib_DiagonalizeSymStack(kineticEvals, devSym, kineticEvecs)
    if (devSym > tolerance) error stop "SysKinetic_SolveShifted (ylm.laplacian.fedvr): T is not symmetric in the grid metric"

    write (*, '(A,I0,A,ES12.4,A,ES12.4,A,ES10.2)') " kinetic blocks for SysKinetic_SolveShifted, l <= ", lmax, &
      ": eigenvalues from", minval(kineticEvals), " to", maxval(kineticEvals), "  asymmetry", devSym

    !$omp atomic write
    solveShiftedReadyQ = .true.

  end subroutine

end submodule
