! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Exchange-correlation functionals on 1D grids (`sysXc.linear`).
!>
!> @details
!> The xc points are the nodes of grid.linear (const or fedvr) with the
!> weights Grid_Linear_weights; densities, potentials and the Hartree source
!> are the nodal defaults of SysXc (S_SysXc). Only the one-dimensional LDA
!> functionals of Libxc qualify (SysXc_dimension = 1; Libxc has no 1D
!> gradient functionals), e.g. for the soft-Coulomb interaction
!> 1/sqrt(x² + 1) of sysInteraction.linear.softYukawa with softening1 = 1:
!>
!> @code{.json}
!> "sysXc": { "linear": { "libxc": {
!>   "functionals": ["lda_x_1d_soft", "lda_c_1d_csc"],
!>   "extParams": { "lda_x_1d_soft": [1.0], "lda_c_1d_csc": [1, 1.0] } } } }
!> @endcode
module M_SysXc_Linear
  use M_Utils_Types

  implicit none

  interface
    !> @brief Bind the 1D xc points, dispatch to the functional.
    module subroutine SysXc_Linear_Fabricate
    end subroutine
  end interface

end module
