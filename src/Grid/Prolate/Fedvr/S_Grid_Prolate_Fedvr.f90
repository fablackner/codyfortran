! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Prolate_Fedvr) S_Grid_Prolate_Fedvr

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Prolate_Fedvr_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Prolate

    call Say_Fabricate("grid.prolate.fedvr")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Grid_Prolate_Fedvr_nElements = Json_Get("grid.prolate.fedvr.nElements", 10)
    Grid_Prolate_Fedvr_nLocals = Json_Get("grid.prolate.fedvr.nLocals", 7)

    ! Total xi points: nE*(nLocals-1) + 1 minus the excluded ximax endpoint
    Grid_Prolate_nXi = Grid_Prolate_Fedvr_nElements * (Grid_Prolate_Fedvr_nLocals - 1)

    ! Update total grid size to include the azimuthal channels
    Grid_nPoints = Grid_Prolate_nXi * Grid_Prolate_nEta * (2 * Grid_Prolate_mmax + 1)

    ! Set the setup procedure
    Grid_Prolate_Setup => Setup

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use stdlib_quadrature, only: gauss_legendre_lobatto
    use M_Utils_Say
    use M_Utils_Fedvr
    use M_Utils_DerivativeFedvr
    use M_Grid
    use M_Grid_Prolate

    integer(I32) :: iE, kLocal, iLocal1, iLocal2, iLocalStart, iLocalEnd, iGlobal1, iGlobal2
    real(R64), allocatable :: elementPoints(:), elementWeights(:)

    call Say_Setup("grid.prolate.fedvr")

    allocate (Grid_Prolate_xiPoints(Grid_Prolate_nXi))
    allocate (Grid_Prolate_xiWeights(Grid_Prolate_nXi))

    ! xi = 1 is kept (sigma orbitals are nonzero on the internuclear axis),
    ! ximax is excluded (Dirichlet box boundary)
    call Fedvr_CreateCtx(Grid_Prolate_Fedvr_fedvrCtx, &
                         1.0_R64, &
                         Grid_Prolate_ximax, &
                         Grid_Prolate_Fedvr_nElements, &
                         Grid_Prolate_Fedvr_nLocals, &
                         .false., &
                         .true.)

    call DerivativeFedvr_CreateCtx(Grid_Prolate_Fedvr_derivativeCtx, Grid_Prolate_Fedvr_fedvrCtx)

    Grid_Prolate_xiPoints(:) = Grid_Prolate_Fedvr_fedvrCtx % points(1:Grid_Prolate_nXi)
    Grid_Prolate_xiWeights(:) = Grid_Prolate_Fedvr_fedvrCtx % weights(1:Grid_Prolate_nXi)

    !------------------------------------
    ! Sturm-Liouville matrix S_xi(i,j) = sum_k w_k (xi_k^2 - 1) D(k,i) D(k,j)
    !------------------------------------

    ! Accumulated element by element in weak form; the quadrature sum runs over
    ! all local nodes (including the excluded ximax endpoint), while the basis
    ! indices skip excluded endpoints. The reference derivative matrix lives on
    ! [-0.5, 0.5], so physical derivatives carry a factor 1/elementSize.
    allocate (Grid_Prolate_xiKinMatrix(Grid_Prolate_nXi, Grid_Prolate_nXi))
    Grid_Prolate_xiKinMatrix = 0.0_R64

    allocate (elementPoints(Grid_Prolate_Fedvr_nLocals))
    allocate (elementWeights(Grid_Prolate_Fedvr_nLocals))

    do iE = 1, Grid_Prolate_Fedvr_nElements
      associate (element => Grid_Prolate_Fedvr_fedvrCtx % elements(iE), &
                 dRef => Grid_Prolate_Fedvr_derivativeCtx % firstOrderMatrix, &
                 nLoc => Grid_Prolate_Fedvr_nLocals)

        call gauss_legendre_lobatto(elementPoints, elementWeights, [element % leftBound, element % rightBound])

        iLocalStart = 1
        iLocalEnd = nLoc
        if (iE .eq. Grid_Prolate_Fedvr_nElements) iLocalEnd = nLoc - 1  ! ximax excluded

        do iLocal2 = iLocalStart, iLocalEnd
          iGlobal2 = element % iStart + (iLocal2 - iLocalStart)
          do iLocal1 = iLocalStart, iLocalEnd
            iGlobal1 = element % iStart + (iLocal1 - iLocalStart)

            do kLocal = 1, nLoc
              Grid_Prolate_xiKinMatrix(iGlobal1, iGlobal2) = Grid_Prolate_xiKinMatrix(iGlobal1, iGlobal2) + &
                elementWeights(kLocal) * (elementPoints(kLocal)**2 - 1.0_R64) * &
                dRef(kLocal, iLocal1) * dRef(kLocal, iLocal2) / element % size**2
            end do

          end do
        end do

      end associate
    end do

    deallocate (elementPoints, elementWeights)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
