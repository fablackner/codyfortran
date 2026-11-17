! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Exchange-correlation functionals on the finite-element grid (`sysXc.fem`).
!>
!> @details
!> The xc points are the interior nodes of grid.fem (SysXc_nPoints =
!> Grid_nPoints) with the lumped-mass weights Grid_weights; densities,
!> potentials and the Hartree source are the nodal defaults of SysXc
!> (S_SysXc). Gradient functionals use the nodal gradient of
!> Grid_Fem_ApplyGradient (element derivatives averaged at shared nodes) and
!> its W-adjoint Grid_Fem_AddGradientAdjoint, so the GGA potential is the
!> exact derivative of the discrete E_xc = Σ_i w_i e(ρ_i, |Gρ|²_i).
!>
!> Requires grid.fem; the Hartree potential uses sysInteraction.fem.coulomb.
!>
!> ## JSON Configuration
!>
!> @code{.json}
!> "sysXc": { "fem": { "libxc": { "functionals": ["gga_x_pbe", "gga_c_pbe"] } } }
!> @endcode
module M_SysXc_Fem
  use M_Utils_Types

  implicit none

  interface
    !> @brief Bind the Fem xc points and gradients, dispatch to the functional.
    module subroutine SysXc_Fem_Fabricate
    end subroutine
  end interface

end module
