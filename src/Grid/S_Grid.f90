! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid) S_Grid

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid_Linear
    use M_Grid_Square
    use M_Grid_Polar
    use M_Grid_Spherical
    use M_Grid_Ylm
    use M_Grid_Lattice

    call Say_Fabricate("grid")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Grid_Orthonormalize => Orthonormalize
    Grid_ProjectOnSubspace => ProjectOnSubspace

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("grid.linear")) then
      call Grid_Linear_Fabricate

    else if (Json_GetExistence("grid.square")) then
      call Grid_Square_Fabricate

    else if (Json_GetExistence("grid.polar")) then
      call Grid_Polar_Fabricate

    else if (Json_GetExistence("grid.spherical")) then
      call Grid_Spherical_Fabricate

    else if (Json_GetExistence("grid.ylm")) then
      call Grid_Ylm_Fabricate

    else if (Json_GetExistence("grid.lattice")) then
      call Grid_Lattice_Fabricate

    else
      error stop "grid is missing one of: linear, square, cubic, polar, spherical, ylm, lattice, fem"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Orthonormalize(fSet)

    complex(R64), intent(inout), contiguous  :: fSet(:, :)

    integer(I32) :: i1, j1, nSet, nG
    real(R64) :: norm
    complex(R64) :: ovlp

    nG = Grid_nPoints
    nSet = size(fSet, 2)

    ! schmidt orthogonalization
    do i1 = 1, nSet
      do j1 = 1, i1 - 1

        ovlp = Grid_InnerProduct(fSet(:, j1), fSet(:, i1))

        fSet(:, i1) = fSet(:, i1) - fSet(:, j1) * ovlp
      end do

      norm = real(Grid_InnerProduct(fSet(:, i1), fSet(:, i1)), kind=R64)
      fSet(:, i1) = fSet(:, i1) / sqrt(norm)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ProjectOnSubspace(fSetInOut, fSet)

    complex(R64), intent(inout), contiguous  :: fSetInOut(:, :)
    complex(R64), intent(in), contiguous     :: fSet(:, :)

    integer(I32) :: i1, j1, nSet, nG
    complex(R64), allocatable :: inp(:, :)

    nG = Grid_nPoints
    nSet = size(fSet, 2)

    allocate (inp(nSet, nSet))

    do j1 = 1, nSet
      do i1 = 1, nSet

        inp(i1, j1) = Grid_InnerProduct(fSet(:, i1), fSetInOut(:, j1))

      end do
    end do

    do j1 = 1, nSet
      do i1 = 1, nSet

        ! ProjectOnSubspace is very delecate and tends to get lost since conservation of orthonormality
        ! requires orthonormality and numerical errors in combination with the nonlinearity of the correlation operator
        ! coupling to the coeffs evolution tends to break this delecate balance.

        fSetInOut(:, j1) = fSetInOut(:, j1) - inp(i1, j1) * fSet(:, i1)

      end do
    end do

  end subroutine

end submodule
