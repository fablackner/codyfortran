! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Lattice_NearestNeighbor.f90
!> @brief Nearest-neighbor tight-binding kinetic operator on 3D lattice.
!>
!> @details
!> Implements the standard tight-binding kinetic energy with direction-dependent
!> hopping amplitudes (hoppX, hoppY, hoppZ). The operator sums hopping
!> contributions from all 6 nearest neighbors (±x, ±y, ±z).
!>
!> Physics
!> -------
!> The tight-binding kinetic operator is:
!>
!>     T̂|i⟩ = −t_x (|i+x⟩ + |i−x⟩) − t_y (|i+y⟩ + |i−y⟩) − t_z (|i+z⟩ + |i−z⟩)
!>
!> where t_x, t_y, t_z are hopping amplitudes. The negative sign convention
!> ensures positive bandwidth. Boundary conditions (periodic or hard-wall) are
!> inherited from the Grid_Lattice module.
!>
!> Connection to Continuum
!> -----------------------
!> For a uniform lattice with spacing a, this approximates the Laplacian:
!>
!>     ∇² ψ ≈ (1/a²) Σ_d [ψ(r+a·d̂) + ψ(r−a·d̂) − 2ψ(r)]
!>
!> The hopping t = 1/(2ma²) in atomic units recovers the correct kinetic energy.
!>
!> @see M_Grid_Lattice for boundary condition settings
module M_SysKinetic_Lattice_NearestNeighbor
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Configure nearest-neighbor hopping and bind the operator.
    !>
    !> JSON keys (with defaults):
    !>   - `sysKinetic.lattice.nearestNeighbor.hoppX` (1.0): x-hopping amplitude
    !>   - `sysKinetic.lattice.nearestNeighbor.hoppY` (1.0): y-hopping amplitude
    !>   - `sysKinetic.lattice.nearestNeighbor.hoppZ` (1.0): z-hopping amplitude
    module subroutine SysKinetic_Lattice_NearestNeighbor_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data (configurable hopping amplitudes)
  !=============================================================================

  !> Hopping amplitude along x-direction [default: 1.0]
  real(R64) :: SysKinetic_Lattice_NearestNeighbor_hoppX = 1.0_R64

  !> Hopping amplitude along y-direction [default: 1.0]
  real(R64) :: SysKinetic_Lattice_NearestNeighbor_hoppY = 1.0_R64

  !> Hopping amplitude along z-direction [default: 1.0]
  real(R64) :: SysKinetic_Lattice_NearestNeighbor_hoppZ = 1.0_R64

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
