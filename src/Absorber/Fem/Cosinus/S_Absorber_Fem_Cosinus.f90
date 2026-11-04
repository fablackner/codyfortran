! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the cosine-profile boundary-layer absorber for the
!> FEM grid.
!>
!> @details
!> The 3D mask is the product of the per-axis layer profiles; it is built
!> once in Setup and applied by element-wise multiplication.
submodule(M_Absorber_Fem_Cosinus) S_Absorber_Fem_Cosinus

  implicit none

  !-----------------------------------------------------------------------------
  ! local data
  !-----------------------------------------------------------------------------

  !> Mask M at every interior node (Grid_nPoints), values in (0, 1].
  real(R64), allocatable :: maskFunction(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Read and validate width and order, bind the procedure pointers.
  module subroutine Absorber_Fem_Cosinus_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Absorber
    use M_Grid_Fem

    real(R64) :: smallestSide

    call Say_Fabricate("absorber.fem.cosinus")

    !------------------------------------
    ! read configuration parameters
    !------------------------------------

    smallestSide = minval(Grid_Fem_axes(:) % upperBound - Grid_Fem_axes(:) % lowerBound)

    Absorber_Fem_Cosinus_width = Json_Get("absorber.fem.cosinus.width", 0.2_R64 * smallestSide)
    Absorber_Fem_Cosinus_order = Json_Get("absorber.fem.cosinus.order", 6)

    if (.not. (Absorber_Fem_Cosinus_width > 0.0_R64 .and. Absorber_Fem_Cosinus_width <= 0.5_R64 * smallestSide)) then
      error stop "absorber.fem.cosinus.width must satisfy 0 < width <= (smallest box side)/2"
    end if

    if (Absorber_Fem_Cosinus_order < 1) error stop "absorber.fem.cosinus.order must be >= 1"

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    Absorber_Setup => Setup
    Absorber_ApplyAbsorber => ApplyAbsorber

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Build the mask M = m_x(x) m_y(y) m_z(z) at the interior nodes.
  subroutine Setup
    !$ use omp_lib, only: omp_in_parallel
    use M_Utils_Say
    use M_Grid
    use M_Grid_Fem

    real(R64), allocatable :: maskX(:), maskY(:), maskZ(:)
    integer(I32) :: ix, iy, iz

    call Say_Setup("absorber.fem.cosinus")

    call FillAxisMask(maskX, Grid_Fem_axes(1))
    call FillAxisMask(maskY, Grid_Fem_axes(2))
    call FillAxisMask(maskZ, Grid_Fem_axes(3))

    if (allocated(maskFunction)) deallocate (maskFunction)
    allocate (maskFunction(Grid_nPoints))

    !$omp parallel do default(shared) private(ix, iy, iz) schedule(static) if (.not. omp_in_parallel())
    do iz = 1, Grid_Fem_nz
      do iy = 1, Grid_Fem_ny
        do ix = 1, Grid_Fem_nx
          maskFunction(Grid_Fem_Index(ix, iy, iz)) = maskX(ix) * maskY(iy) * maskZ(iz)
        end do
      end do
    end do
    !$omp end parallel do

    write (*, '(A, F10.4, A, I0, A, I0, A, I0)') "  absorber.fem.cosinus: width ", Absorber_Fem_Cosinus_width, &
      ", layer nodes per axis (x, y, z): ", count(maskX < 1.0_R64), ", ", count(maskY < 1.0_R64), ", ", &
      count(maskZ < 1.0_R64)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Layer profile of one axis at its interior nodes:
  !> cos^(1/order)(pi/2 (width - d)/width) for d < width, 1 otherwise, d the
  !> distance to the nearer face.
  subroutine FillAxisMask(mask, axis)
    use M_Utils_Constants, only: PI
    use M_Grid_Fem

    real(R64), allocatable, intent(out) :: mask(:)
    type(T_Grid_Fem_Axis), intent(in)   :: axis

    real(R64) :: width, power, distance, profile
    integer(I32) :: i

    width = Absorber_Fem_Cosinus_width
    power = 1.0_R64 / real(Absorber_Fem_Cosinus_order, R64)

    allocate (mask(axis % nPoints))

    do i = 1, axis % nPoints
      distance = min(axis % points(i) - axis % lowerBound, axis % upperBound - axis % points(i))

      if (distance < width) then
        ! the faces (distance 0) are Dirichlet nodes, so the profile is > 0
        profile = max(cos(0.5_R64 * PI * (width - distance) / width), 0.0_R64)
        mask(i) = profile**power
      else
        mask(i) = 1.0_R64
      end if
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Multiply every orbital column by the mask, in place.
  !>
  !> @param[in,out] orbs  Orbitals (Grid_nPoints, nOrbitals).
  subroutine ApplyAbsorber(orbs)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(inout), contiguous :: orbs(:, :)

    integer(I32) :: i, iOrb

    if (.not. allocated(maskFunction)) then
      error stop "absorber.fem.cosinus: call Absorber_Setup before Absorber_ApplyAbsorber"
    end if
    if (size(orbs, 1) .ne. size(maskFunction)) then
      error stop "absorber.fem.cosinus: orbitals must have Grid_nPoints rows"
    end if

    do iOrb = 1, size(orbs, 2)
      !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
      do i = 1, size(maskFunction)
        orbs(i, iOrb) = orbs(i, iOrb) * maskFunction(i)
      end do
      !$omp end parallel do
    end do

  end subroutine

end submodule
