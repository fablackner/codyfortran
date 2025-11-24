! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Propagator_EigenExpansion) S_Propagator_EigenExpansion

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Propagator_EigenExpansion_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Propagator

    call Say_Fabricate("propagator.eigenExpansion")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Propagator_Propagate => Propagate
    Propagator_Setup => Setup

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup()
    use M_Utils_Say
    use M_DiagonalizerList

    call Say_Setup("propagator.eigenExpansion")

    call DiagonalizerList(1) % e % Diagonalize(0.0_R64, .true.)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Propagate(state, t0, t1)
    use M_Utils_Constants
    use M_DiagonalizerList

    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in)    :: t0      ! Starting time
    real(R64), intent(in)    :: t1      ! Target ending time

    integer(I32) :: i, nEvals
    real(R64) :: dt
    complex(R64), allocatable :: c(:)

    dt = t1 - t0
    nEvals = DiagonalizerList(1) % e % nFound

    allocate (c(nEvals))

    do i = 1, nEvals
      c(i) = dot_product(DiagonalizerList(1) % e % evecs(:, i), state)
    end do

    state(:) = 0.0_R64
    do i = 1, nEvals
      associate (diagonalizer => DiagonalizerList(1) % e)
        state(:) = state(:) + exp(-IU * diagonalizer % evals(i) * dt) * c(i) * diagonalizer % evecs(:, i)
      end associate
    end do

  end subroutine

end submodule
