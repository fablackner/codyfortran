! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for M_OrbsInit_Fem.
!>
!> @details
!> Fabrication and the generic InitializeOrb, which samples the model's
!> InitFunction at the nodes and normalizes.
submodule(M_OrbsInit_Fem) S_OrbsInit_Fem

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Check the grid, bind InitializeOrb and dispatch to the model.
  module subroutine OrbsInit_Fem_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit
    use M_OrbsInit_Fem_Sto

    call Say_Fabricate("orbsInit.fem")

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.fem")) then
      error stop "grid.fem is required for orbsInit.fem"
    end if

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_InitializeOrb => InitializeOrb

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("orbsInit.fem.sto")) then
      call OrbsInit_Fem_Sto_Fabricate

    else
      error stop "orbsInit.fem is missing one of: sto"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Initialize one orbital on the FEM nodes.
!>
!> @details
!> Samples OrbsInit_Fem_InitFunction at every interior node (OpenMP over the
!> z planes outside parallel regions) and normalizes with Grid_InnerProduct.
!> The orbital is not orthogonalized against others.
!>
!> @param[out] orb  Orbital on the grid (Grid_nPoints), normalized.
!> @param[in]  ind  Orbital index (1..OrbsInit_Fem_nOrbitals).
!> @param[in]  bt_  Optional body type (absent for method.sb), passed through.
  subroutine InitializeOrb(orb, ind, bt_)
    !$ use omp_lib, only: omp_in_parallel
    use M_Grid
    use M_Grid_Fem

    complex(R64), intent(out), contiguous :: orb(:)
    integer(I32), intent(in)              :: ind
    integer(I32), intent(in), optional    :: bt_

    character(len=200) :: message
    real(R64) :: x, y, z, norm
    integer(I32) :: ix, iy, iz

    if (ind < 1 .or. ind > OrbsInit_Fem_nOrbitals) then
      write (message, '(A, I0, A, I0, A)') "orbsInit.fem: orbital ", ind, " requested, but the model defines ", &
        OrbsInit_Fem_nOrbitals, " orbitals"
      error stop trim(message)
    end if

    if (size(orb) .ne. Grid_nPoints) error stop "orbsInit.fem: the orbital must have Grid_nPoints entries"

    ! Sample the initialization function at every node
    !$omp parallel do default(shared) private(ix, iy, iz, x, y, z) schedule(static) if (.not. omp_in_parallel())
    do iz = 1, Grid_Fem_nz
      z = Grid_Fem_axes(3) % points(iz)
      do iy = 1, Grid_Fem_ny
        y = Grid_Fem_axes(2) % points(iy)
        do ix = 1, Grid_Fem_nx
          x = Grid_Fem_axes(1) % points(ix)
          orb(Grid_Fem_Index(ix, iy, iz)) = cmplx(OrbsInit_Fem_InitFunction(x, y, z, ind, bt_=bt_), 0.0_R64, kind=R64)
        end do
      end do
    end do
    !$omp end parallel do

    ! Normalize with the grid metric
    norm = real(Grid_InnerProduct(orb, orb), kind=R64)

    if (.not. (norm > 0.0_R64)) then
      write (message, '(A, I0, A)') "orbsInit.fem: orbital ", ind, " vanishes on the grid"
      error stop trim(message)
    end if

    orb(:) = orb(:) / sqrt(norm)

  end subroutine

end submodule
