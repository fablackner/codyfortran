! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default numerical implementation for lattice harmonic potential.
!>
!> Evaluates a separable 3D harmonic trap on a discrete lattice:
!>   V(i,j,k) = ½ωₓ(i-x₀)² + ½ωᵧ(j-y₀)² + ½ω_z(k-z₀)²
!>
!> The trap center (x₀,y₀,z₀) and frequencies (ωₓ,ωᵧ,ω_z) are read from the
!> parent module `M_SysPotential_Lattice_Harmonic`.
submodule(M_SysPotential_Lattice_Harmonic_StdImpl) S_SysPotential_Lattice_Harmonic_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Lattice_Harmonic_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Lattice

    call Say_Fabricate("sysPotential.lattice.harmonic.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_FillExternalPotential => FillExternalPotential

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Fill the external potential array with harmonic trap values.
  !>
  !> Iterates over all lattice sites and computes V = ½ω²(r-r₀)² for each axis.
  !> The time and body-type arguments are unused (time-independent, body-type-independent).
  subroutine FillExternalPotential(externalPotential, time, bt_)
    use M_Utils_UnusedVariables
    use M_SysPotential_Lattice_Harmonic
    use M_Grid
    use M_Grid_Lattice

    complex(R64), intent(out), allocatable :: externalPotential(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: ix, iy, iz, iGrid
    real(R64) :: omegaX, omegaY, omegaZ
    real(R64) :: positionX, positionY, positionZ

    if (.false.) call UnusedVariables_Mark(time, bt_)

    if (.not. allocated(externalPotential)) allocate (externalPotential(Grid_nPoints))

    omegaX = SysPotential_Lattice_Harmonic_omegaX
    omegaY = SysPotential_Lattice_Harmonic_omegaY
    omegaZ = SysPotential_Lattice_Harmonic_omegaZ
    positionX = SysPotential_Lattice_Harmonic_positionX
    positionY = SysPotential_Lattice_Harmonic_positionY
    positionZ = SysPotential_Lattice_Harmonic_positionZ

    do iz = 1, Grid_Lattice_zSize
      do iy = 1, Grid_Lattice_ySize
        do ix = 1, Grid_Lattice_xSize
          iGrid = Grid_Lattice_code(ix, iy, iz)
          externalPotential(iGrid) = 0.5_R64 * ( &
                                     omegaX * (ix - positionX)**2 + &
                                     omegaY * (iy - positionY)**2 + &
                                     omegaZ * (iz - positionZ)**2)
        end do
      end do
    end do

  end subroutine

end submodule
