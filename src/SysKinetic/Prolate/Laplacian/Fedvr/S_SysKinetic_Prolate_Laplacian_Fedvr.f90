! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Prolate_Laplacian_Fedvr.f90
!> @brief Implementation of the FEDVR kinetic operator for prolate channels.
submodule(M_SysKinetic_Prolate_Laplacian_Fedvr) S_SysKinetic_Prolate_Laplacian_Fedvr

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind FEDVR channel operator, validate grid requirements.
  module subroutine SysKinetic_Prolate_Laplacian_Fedvr_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Prolate

    call Say_Fabricate("sysKinetic.prolate.laplacian.fedvr")

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    SysKinetic_Prolate_MultiplyWithChannelKineticOp => MultiplyWithChannelKineticOp

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.prolate.fedvr")) then
      error stop "grid.prolate.fedvr is required for sysKinetic.prolate.laplacian.fedvr"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply the kinetic operator to one azimuthal channel.
  !>
  !> @details
  !> Algorithm on the (nXi, nEta) channel matrix f:
  !>   1. Sturm-Liouville terms via the grid's symmetric DVR matrices:
  !>      [d/dxi (xi^2-1) d/dxi f]  = -(1/w_xi)  S_xi  f    (rows)
  !>      [d/deta (1-eta^2) d/deta f] = -(1/w_eta) f S_eta  (columns)
  !>   2. Centrifugal term m^2 (1/(xi^2-1) + 1/(1-eta^2)) f for m /= 0
  !>   3. Common prefactor -1/(2 M a^2 (xi^2 - eta^2))
  !>
  !> For m /= 0 the channel physically vanishes on the internuclear axis like
  !> (xi^2-1)^(|m|/2); the xi = 1 grid point (where the centrifugal term
  !> diverges) is therefore treated as a per-channel Dirichlet point: its
  !> input and output values are zeroed.
  subroutine MultiplyWithChannelKineticOp(dOrbM, orbM, m, time, bt_)
    use M_Utils_UnusedVariables
    use M_Grid
    use M_Grid_Prolate
    use M_SysKinetic_Prolate_Laplacian

    complex(R64), intent(out) :: dOrbM(:)
    complex(R64), intent(in)  :: orbM(:)
    integer(I32), intent(in)  :: m
    real(R64), intent(in)     :: time
    integer(I32), intent(in), optional :: bt_

    complex(R64), allocatable :: f(:, :), res(:, :)
    integer(I32) :: bt, iXi, jEta
    real(R64) :: centrifugal

    if (.false.) call UnusedVariables_Mark(time)

    bt = 1
    if (present(bt_)) bt = bt_

    allocate (f(Grid_Prolate_nXi, Grid_Prolate_nEta))
    allocate (res(Grid_Prolate_nXi, Grid_Prolate_nEta))

    f = reshape(orbM, [Grid_Prolate_nXi, Grid_Prolate_nEta])

    ! Dirichlet at xi = 1 for m /= 0 (channel vanishes on the axis)
    if (m .ne. 0) f(1, :) = (0.0_R64, 0.0_R64)

    ! Sturm-Liouville derivative terms (positive semi-definite convention:
    ! -d/dxi (xi^2-1) d/dxi corresponds to +S_xi/w_xi)
    res = matmul(Grid_Prolate_xiKinMatrix, f)
    do iXi = 1, Grid_Prolate_nXi
      res(iXi, :) = res(iXi, :) / Grid_Prolate_xiWeights(iXi)
    end do

    f = matmul(f, Grid_Prolate_etaKinMatrix) ! S_eta is symmetric
    do jEta = 1, Grid_Prolate_nEta
      res(:, jEta) = res(:, jEta) + f(:, jEta) / Grid_Prolate_etaWeights(jEta)
    end do

    ! Restore f for the centrifugal term
    f = reshape(orbM, [Grid_Prolate_nXi, Grid_Prolate_nEta])

    do jEta = 1, Grid_Prolate_nEta
      do iXi = 1, Grid_Prolate_nXi

        if (m .ne. 0 .and. iXi > 1) then
          centrifugal = m * m * (1.0_R64 / (Grid_Prolate_xiPoints(iXi)**2 - 1.0_R64) + &
                                 1.0_R64 / (1.0_R64 - Grid_Prolate_etaPoints(jEta)**2))
          res(iXi, jEta) = res(iXi, jEta) + centrifugal * f(iXi, jEta)
        end if

        ! Common prefactor: T = +1/(2 M a^2 (xi^2-eta^2)) [S terms + centrifugal]
        res(iXi, jEta) = res(iXi, jEta) / &
                         (2.0_R64 * SysKinetic_Prolate_Laplacian_bodyMass(bt) * Grid_Prolate_a**2 * &
                          (Grid_Prolate_xiPoints(iXi)**2 - Grid_Prolate_etaPoints(jEta)**2))
      end do
    end do

    if (m .ne. 0) res(1, :) = (0.0_R64, 0.0_R64)

    dOrbM(:) = reshape(res, [Grid_Prolate_nSpatial])

    deallocate (f, res)

  end subroutine

end submodule
