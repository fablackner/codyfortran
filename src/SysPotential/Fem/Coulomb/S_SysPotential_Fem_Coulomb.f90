! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Implementation submodule for the static nuclear Coulomb attraction on
!> grid.fem.
submodule(M_SysPotential_Fem_Coulomb) S_SysPotential_Fem_Coulomb

  implicit none

  !-----------------------------------------------------------------------------
  ! local data
  !-----------------------------------------------------------------------------

  !> Static (time- and body-type-independent) potential at the nodes
  !> (Grid_nPoints), built once in Setup and copied by FillExternalPotential.
  complex(R64), allocatable :: staticPotential(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Fem_Coulomb_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid_Fem
    use M_SysPotential, only: SysPotential_Setup, SysPotential_FillExternalPotential, &
                              SysPotential_timeIndependentQ, SysPotential_bodyTypeIndependentQ

    character(len=200) :: message
    integer(I32) :: iCenter

    call Say_Fabricate("sysPotential.fem.coulomb")

    !------------------------------------
    ! charges and smearing radii, one per grid.fem center
    !------------------------------------

    if (.not. Json_GetExistence("sysPotential.fem.coulomb.charges")) then
      error stop "sysPotential.fem.coulomb.charges is required (one nuclear charge per grid.fem center, 0 for a ghost)"
    end if

    SysPotential_Fem_Coulomb_charges = Json_Get("sysPotential.fem.coulomb.charges", [real(R64) ::])
    SysPotential_Fem_Coulomb_smearing = Json_Get("sysPotential.fem.coulomb.smearing", &
                                                 [(0.0_R64, iCenter=1, Grid_Fem_nCenters)])

    if (size(SysPotential_Fem_Coulomb_charges) .ne. Grid_Fem_nCenters) then
      write (message, '(A, I0, A, I0, A)') "sysPotential.fem.coulomb.charges has ", size(SysPotential_Fem_Coulomb_charges), &
        " entries, grid.fem has ", Grid_Fem_nCenters, " centers (one charge per center)"
      error stop trim(message)
    end if

    if (size(SysPotential_Fem_Coulomb_smearing) .ne. Grid_Fem_nCenters) then
      write (message, '(A, I0, A, I0, A)') "sysPotential.fem.coulomb.smearing has ", size(SysPotential_Fem_Coulomb_smearing), &
        " entries, grid.fem has ", Grid_Fem_nCenters, " centers (one radius per center)"
      error stop trim(message)
    end if

    if (any(.not. (SysPotential_Fem_Coulomb_smearing >= 0.0_R64))) then
      error stop "sysPotential.fem.coulomb.smearing must be >= 0 (0: Galerkin point nucleus)"
    end if

    if (any(.not. (abs(SysPotential_Fem_Coulomb_charges) < huge(1.0_R64)))) then
      error stop "sysPotential.fem.coulomb.charges must be finite"
    end if

    SysPotential_Fem_Coulomb_nuclearRepulsion = NuclearRepulsion(SysPotential_Fem_Coulomb_charges, Grid_Fem_centers)

    !------------------------------------
    ! flags and procedure pointers
    !------------------------------------

    SysPotential_timeIndependentQ = .true.
    SysPotential_bodyTypeIndependentQ = .true.

    SysPotential_Setup => Setup
    SysPotential_FillExternalPotential => FillExternalPotential

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure module function SysPotential_Fem_Coulomb_SmearedPotential(charge, smearing, r) result(res)
    use M_Utils_Constants, only: PI

    real(R64), intent(in) :: charge
    real(R64), intent(in) :: smearing
    real(R64), intent(in) :: r
    real(R64) :: res

    ! erf(x)/x = 2/sqrt(pi) (1 - x^2/3 + ...): the limit is exact to 1e-24
    ! relative below x = 1e-12
    if (r < 1.0e-12_R64 * smearing) then
      res = -2.0_R64 * charge / (sqrt(PI) * smearing)
    else
      res = -charge * erf(r / smearing) / r
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Build the static potential once: all point nuclei in one Galerkin solve
  !> (charges -Z), then the nodal values of the smeared nuclei; ghost centers
  !> are skipped. Prints the nuclear repulsion and the potential minimum.
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Grid_Fem

    real(R64), allocatable :: pointCharges(:), pointPositions(:, :)
    integer(I32) :: iCenter, nPoint, nSmeared, nGhost, iPoint

    call Say_Setup("sysPotential.fem.coulomb")

    if (allocated(staticPotential)) deallocate (staticPotential)
    allocate (staticPotential(Grid_nPoints))

    nPoint = 0
    nSmeared = 0
    nGhost = 0
    do iCenter = 1, Grid_Fem_nCenters
      if (GhostQ(iCenter)) then
        nGhost = nGhost + 1
      else if (SysPotential_Fem_Coulomb_smearing(iCenter) > 0.0_R64) then
        nSmeared = nSmeared + 1
      else
        nPoint = nPoint + 1
      end if
    end do

    ! Galerkin point nuclei: one solve for all of them
    if (nPoint > 0) then
      allocate (pointCharges(nPoint), pointPositions(3, nPoint))
      iPoint = 0
      do iCenter = 1, Grid_Fem_nCenters
        if (GhostQ(iCenter) .or. SysPotential_Fem_Coulomb_smearing(iCenter) > 0.0_R64) cycle
        iPoint = iPoint + 1
        pointCharges(iPoint) = -SysPotential_Fem_Coulomb_charges(iCenter)
        pointPositions(:, iPoint) = Grid_Fem_centers(:, iCenter)
      end do

      call Grid_Fem_SolvePoissonPointCharges(staticPotential, pointCharges, pointPositions)
    else
      staticPotential(:) = (0.0_R64, 0.0_R64)
    end if

    ! Smeared nuclei: nodal values
    do iCenter = 1, Grid_Fem_nCenters
      if (GhostQ(iCenter) .or. .not. (SysPotential_Fem_Coulomb_smearing(iCenter) > 0.0_R64)) cycle
      call AddSmearedNucleus(staticPotential, SysPotential_Fem_Coulomb_charges(iCenter), &
                             SysPotential_Fem_Coulomb_smearing(iCenter), Grid_Fem_centers(:, iCenter), &
                             Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz)
    end do

    write (*, '(A, I0, A, I0, A, I0, A)') "  sysPotential.fem.coulomb: ", nPoint, " point, ", nSmeared, &
      " smeared and ", nGhost, " ghost centers"
    write (*, '(A, F22.12)') "  nuclear repulsion                  = ", SysPotential_Fem_Coulomb_nuclearRepulsion
    write (*, '(A, F22.12)') "  minimum of the static potential    = ", &
      minval(real(staticPotential, kind=R64))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Copy the static potential built in Setup (time and bt_ are irrelevant).
  subroutine FillExternalPotential(externalPotential, time, bt_)
    use M_Utils_UnusedVariables

    complex(R64), intent(out), allocatable :: externalPotential(:)
    real(R64), intent(in)                  :: time
    integer(I32), intent(in), optional     :: bt_

    if (.false.) call UnusedVariables_Mark(time, bt_)

    if (.not. allocated(staticPotential)) then
      error stop "sysPotential.fem.coulomb: call SysPotential_Setup before SysPotential_FillExternalPotential"
    end if

    externalPotential = staticPotential

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> pot(ix, iy, iz) += -Z erf(r/rc)/r, r the distance of the node from the
  !> nucleus at center (limit -2Z/(sqrt(pi) rc) on the nucleus node).
  subroutine AddSmearedNucleus(pot, charge, smearing, center, nx, ny, nz)
    !$ use omp_lib, only: omp_in_parallel
    use M_Grid_Fem

    integer(I32), intent(in)    :: nx
    integer(I32), intent(in)    :: ny
    integer(I32), intent(in)    :: nz
    complex(R64), intent(inout) :: pot(nx, ny, nz)
    real(R64), intent(in)       :: charge
    real(R64), intent(in)       :: smearing
    real(R64), intent(in)       :: center(3)

    real(R64) :: dx, dy, dz
    integer(I32) :: ix, iy, iz

    !$omp parallel do default(shared) private(ix, iy, iz, dx, dy, dz) schedule(static) if (.not. omp_in_parallel())
    do iz = 1, nz
      dz = Grid_Fem_axes(3) % points(iz) - center(3)
      do iy = 1, ny
        dy = Grid_Fem_axes(2) % points(iy) - center(2)
        do ix = 1, nx
          dx = Grid_Fem_axes(1) % points(ix) - center(1)
          pot(ix, iy, iz) = pot(ix, iy, iz) + &
                            SysPotential_Fem_Coulomb_SmearedPotential(charge, smearing, sqrt(dx * dx + dy * dy + dz * dz))
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Nuclear repulsion sum_(I<J) Z_I Z_J / |R_I - R_J| (pairs with a ghost
  !> center skipped; two charged centers on the same position stop).
  function NuclearRepulsion(charges, centers) result(res)

    real(R64), intent(in) :: charges(:)
    real(R64), intent(in) :: centers(:, :)
    real(R64) :: res

    real(R64) :: distance
    integer(I32) :: i, j

    res = 0.0_R64
    do i = 1, size(charges)
      do j = i + 1, size(charges)
        if (abs(charges(i) * charges(j)) <= 0.0_R64) cycle

        distance = norm2(centers(:, i) - centers(:, j))
        if (.not. (distance > 0.0_R64)) then
          error stop "sysPotential.fem.coulomb: two charged centers coincide (infinite nuclear repulsion)"
        end if

        res = res + charges(i) * charges(j) / distance
      end do
    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Whether center iCenter is a ghost (charge 0).
  function GhostQ(iCenter) result(res)

    integer(I32), intent(in) :: iCenter
    logical :: res

    res = abs(SysPotential_Fem_Coulomb_charges(iCenter)) <= 0.0_R64

  end function

end submodule
