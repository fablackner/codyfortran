! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Simple numerical integration helpers.
module M_Utils_Integrate
  use M_Utils_Types
  implicit none

  private I_Fun

  abstract interface
    function I_Fun(x) result(res)
      import :: R64
      real(R64) :: res
      real(R64) :: x
    end function
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Integrate a scalar function on [xMin, xMax] using a composite Simpson-like rule.
!>
!> The integrator evaluates the function at endpoints, interior grid points,
!> and midpoints between them with weights (1, 2, 4) respectively, and
!> multiplies by dx/6. This is effective for smooth functions on uniform grids.
!>
!> Parameters
!> - fun   [in]: procedure that evaluates f(x).
!> - xMin  [in]: lower bound of integration interval.
!> - xMax  [in]: upper bound of integration interval.
!> - nSteps[in]: number of subintervals; must be >= 1.
!>
!> Returns the approximate integral ∫ f(x) dx over [xMin, xMax].
  function Integrate_RiemannSum(fun, xMin, xMax, nSteps) result(res)
    procedure(I_Fun) :: fun
    real(R64), intent(in) :: xMin, xMax
    integer(I32), intent(in) :: nSteps
    real(R64) :: res
    real(R64) :: dx
    integer(I32) :: ix

    dx = (xMax - xMin) / nSteps
    res = fun(xMin) + fun(xMax)

    do ix = 1, nSteps - 1
      res = res + 2.0_R64 * fun(xMin + ix * dx)
    end do
    do ix = 1, nSteps
      res = res + 4.0_R64 * fun(xMin + (ix - 0.5_R64) * dx)
    end do
    res = res * dx / 6.0_R64
  end function

end module
