! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Prolate.f90
!> @brief Kinetic operators in prolate-spheroidal (two-center) coordinates.
!>
!> @details
!> Fields are stored in azimuthal channels (see M_Grid_Prolate):
!>
!>     psi(xi,eta,phi) = sum_m f_m(xi,eta) e^(i m phi)/sqrt(2 pi)
!>
!> The kinetic operator acts channel-diagonally. With a = R/2 and mass M:
!>
!>     T f_m = -1/(2 M a^2 (xi^2-eta^2)) [ d/dxi (xi^2-1) d/dxi
!>                                        + d/deta (1-eta^2) d/deta
!>                                        - m^2 ( 1/(xi^2-1) + 1/(1-eta^2) ) ] f_m
!>
!> The Sturm-Liouville derivative terms are applied via the symmetric DVR
!> matrices provided by the grid (Grid_Prolate_xiKinMatrix / etaKinMatrix), so
!> T is Hermitian with respect to the metric-weighted inner product.
!>
!> @see M_SysKinetic_Prolate_Laplacian, M_Grid_Prolate
module M_SysKinetic_Prolate
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Fabricate the prolate kinetic backend and bind the channel operator.
    module subroutine SysKinetic_Prolate_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to the kinetic operator for a single azimuthal channel.
  procedure(I_SysKinetic_Prolate_MultiplyWithChannelKineticOp), pointer :: SysKinetic_Prolate_MultiplyWithChannelKineticOp
  abstract interface
    !> @brief Apply the kinetic operator to a single m channel.
    !>
    !> @param[out] dOrbM  Output spatial function (nXi*nEta, xi fastest)
    !> @param[in]  orbM   Input spatial function f_m(xi,eta)
    !> @param[in]  m      Azimuthal quantum number
    !> @param[in]  time   Current simulation time (reserved for extensions)
    !> @param[in]  bt_    Optional body-type index for mass scaling (default: 1)
    subroutine I_SysKinetic_Prolate_MultiplyWithChannelKineticOp(dOrbM, orbM, m, time, bt_)
      import :: I32, R64
      complex(R64), intent(out) :: dOrbM(:)
      complex(R64), intent(in) :: orbM(:)
      integer(I32), intent(in)  :: m
      real(R64), intent(in) :: time
      integer(I32), intent(in), optional :: bt_
    end subroutine
  end interface

end module
