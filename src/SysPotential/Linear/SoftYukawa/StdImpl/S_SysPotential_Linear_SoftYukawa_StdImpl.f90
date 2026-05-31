! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default numerical implementation for linear soft-Yukawa potential.
!>
!> Evaluates a sum of softened Yukawa (screened Coulomb) potentials:
!>   V(x) = -Σₙ qₙ exp(-αₙ|x-xₙ|) / (√((x-xₙ)² + s₁ₙ²) + s₂ₙ)
!>
!> The softening parameters regularize the Coulomb singularity:
!> - s₁ (softening1): Enters under the square root, smooths the 1/r behavior
!> - s₂ (softening2): Additive offset, provides additional short-range control
!> - α (dampening): Yukawa screening length; α=0 gives pure soft-Coulomb
!>
!> This form is commonly used for 1D model atoms and molecules where the
!> true 3D Coulomb potential would be inappropriate.
submodule(M_SysPotential_Linear_SoftYukawa_StdImpl) S_SysPotential_Linear_SoftYukawa_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Linear_SoftYukawa_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Linear

    call Say_Fabricate("sysPotential.linear.softYukawa.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_FillExternalPotential => FillExternalPotential

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Fill the external potential array with soft-Yukawa values.
  !>
  !> Sums contributions from all charge centers. Each center contributes a
  !> softened Coulomb-like potential with optional exponential screening.
  subroutine FillExternalPotential(externalPotential, time, bt_)
    use M_Utils_UnusedVariables
    use M_SysPotential_Linear_SoftYukawa
    use M_Grid
    use M_Grid_Linear

    complex(R64), intent(out), allocatable :: externalPotential(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: i
    real(R64), allocatable :: position(:)
    real(R64), allocatable :: charge(:)
    real(R64), allocatable :: softening1(:)
    real(R64), allocatable :: softening2(:)
    real(R64), allocatable :: dampening(:)
    real(R64), allocatable :: distances(:)

    if (.false.) call UnusedVariables_Mark(time, bt_)

    position = SysPotential_Linear_SoftYukawa_position
    charge = SysPotential_Linear_SoftYukawa_charge
    softening1 = SysPotential_Linear_SoftYukawa_softening1
    softening2 = SysPotential_Linear_SoftYukawa_softening2
    dampening = SysPotential_Linear_SoftYukawa_dampening

    if (.not. allocated(externalPotential)) allocate (externalPotential(Grid_nPoints))
    allocate (distances(size(externalPotential)))

    externalPotential(:) = 0.0_R64

    do i = 1, size(position)
      distances = abs(Grid_Linear_xCoord(:) - position(i))
      externalPotential(:) = externalPotential(:) - charge(i) * exp(-distances * dampening(i)) / &
                             (sqrt(distances**2 + softening1(i)**2) + softening2(i))
    end do

    deallocate (distances)

  end subroutine

end submodule
