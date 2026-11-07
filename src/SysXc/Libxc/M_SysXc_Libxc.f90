! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Exchange-correlation functional of Libxc (`sysXc.<family>.libxc`).
!>
!> @details
!> A linear combination Σ_i c_i f_i of Libxc functionals (LDA, GGA and their
!> global hybrids), evaluated on the xc points of any grid family: the
!> family dispatches here with its JSON path, so one functional
!> implementation serves all grids. Binds `SysXc_EvaluateFunctional` and sets
!> `SysXc_gradientQ` (any GGA) and `SysXc_exactExchangeFraction`
!> (Σ_i c_i a_i of the hybrids, or the JSON override).
!>
!> ## JSON Configuration (under `sysXc.<family>.libxc`)
!>
!> | Key | Type | Default | Description |
!> |-----|------|---------|-------------|
!> | `functionals` | string(:) | required | Libxc names, e.g. ["lda_x", "lda_c_pw"], ["gga_x_pbe", "gga_c_pbe"], ["hyb_gga_xc_pbeh"] (case-insensitive, optional prefix xc_); may be empty (pure exact exchange) |
!> | `coefficients` | real(:) | 1 each | Weights c_i of the functionals |
!> | `densityThreshold` | real | Libxc's | Densities below it count as zero (all functionals) |
!> | `extParams.<name>` | real(:) | Libxc's | External parameters of functional <name> (e.g. lda_x_1d_soft: [beta]) |
!> | `exactExchangeFraction` | real | Σ_i c_i a_i | Fraction of exact exchange (overrides the hybrids') |
!>
!> Refused (error stop): meta-GGAs, range-separated hybrids (CAM, LC), VV10
!> nonlocal correlation, kinetic functionals, functionals without energy or
!> potential, and functionals of another spatial dimension than the grid
!> (`SysXc_dimension`).
!>
!> ## Evaluation
!>
!> In chunks of `SysXc_Libxc_CHUNK_SIZE` points, OpenMP-parallel over the
!> chunks (Libxc's evaluation only reads the functional handles), with
!> per-thread buffers in Libxc's interleaved spin layout.
module M_SysXc_Libxc
  use M_Utils_Types

  implicit none

  !> Points per Libxc call
  integer(I32), parameter :: SysXc_Libxc_CHUNK_SIZE = 2048

  interface
    !> @brief Read the functional below path and bind SysXc_EvaluateFunctional.
    !> @param[in] path  JSON path of the libxc block, e.g. "sysXc.fem.libxc"
    module subroutine SysXc_Libxc_Fabricate(path)
      character(len=*), intent(in) :: path
    end subroutine
  end interface

end module
