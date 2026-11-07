! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Exchange-correlation functionals on Ylm grids (`sysXc.ylm`).
!>
!> @details
!> ## Xc points
!>
!> Every radial point of grid.ylm carries an angular product quadrature:
!> `SysXc_Ylm_nTheta` Gauss-Legendre nodes x_q = cos θ_q times
!> `SysXc_Ylm_nPhi` uniform angles φ_k = 2π(k - 1)/nPhi, laid out as
!> p = iRad + nRad (q - 1) + nRad nTheta (k - 1) with the weights
!> w_p = w_r(iRad) 2π w_q / nPhi (w_r = Grid_Ylm_radialWeights, r² included).
!> An orbital f = Σ_lm f_lm(r) Θ_lm(θ) e^{imφ} (Condon-Shortley Θ_lm, own
!> tables up to the largest l and |m| needed) is evaluated there by one real
!> GEMM per m (channels → nodes) and a sum over m per φ_k; a local potential
!> v acts as (v f)_lm = Σ_p (w_p / w_r) Θ_lm e^{-imφ} v f (nodes → channels),
!> Hermitian in the Ylm metric. Densities and their interaction source,
!> the Hartree potential and the GGA gradients live on the same points.
!>
!> ## Axial mode (nPhi = 1)
!>
!> For orbitals of definite m (method.mb.orbBased.mOfOrbs) or grid.ylm.mmax =
!> 0 the density is axially symmetric: every m block of an orbital is
!> transformed separately, the density is Σ_m |f_m(r, θ)|² (the φ-average,
!> exact for these states) and v(r, θ) multiplies each m block without m
!> coupling. Selected by default in these cases; for other states it is
!> the axial approximation (φ-averaged density, a note is printed).
!>
!> ## General mode (nPhi > 1)
!>
!> Orbitals are synthesized at every φ_k from all their m. Densities have
!> |m| ≤ 2 mmax, so nPhi ≥ 4 mmax + 1 (the default) resolves them and the
!> Hartree potential exactly; v_xc is not band-limited and converges with
!> nPhi and nTheta.
!>
!> ## Hartree potential
!>
!> The interaction source is the channel projection of the density onto the
!> layout of sysInteraction.ylm (l ≤ SysInteraction_Ylm_lmax, |m| ≤
!> SysInteraction_Ylm_mmax, radial metric weights) from the xc points, and
!> the potential channels are evaluated back on them: one Poisson solve per
!> channel, no Gaunt products (exact for nTheta ≥ 2 lmax + 1).
!>
!> ## Gradients (GGA, grid.ylm.fedvr)
!>
!> Components in the spherical frame: ∂_r f = (D(r f) - f)/r with the banded
!> nodal FEDVR derivative D (element derivatives averaged at shared nodes;
!> r f vanishes at the excluded node r = 0), (1/r) ∂_θ f by projection onto
!> the channels l ≤ 2 lmax of a density and differentiation of Θ_lm, and
!> (1/(r sin θ)) ∂_φ f by the discrete Fourier series in φ; their adjoints in
!> the metric of the xc points follow by transposition.
!>
!> ## JSON Configuration
!>
!> | Key | Type | Default | Description |
!> |-----|------|---------|-------------|
!> | `sysXc.ylm.nTheta` | int | 2 lmax + 3 | Gauss-Legendre nodes in cos θ (≥ lmax + 1) |
!> | `sysXc.ylm.nPhi` | int | 1 if mmax = 0 or mOfOrbs, else 4 mmax + 1 | φ points; 1 selects the axial mode |
!>
!> @code{.json}
!> "sysXc": { "ylm": { "libxc": { "functionals": ["lda_x", "lda_c_vwn"] } } }
!> @endcode
module M_SysXc_Ylm
  use M_Utils_Types

  implicit none

  interface
    !> @brief Bind the Ylm xc points and maps, dispatch to the functional.
    module subroutine SysXc_Ylm_Fabricate
    end subroutine
  end interface

  !> Gauss-Legendre nodes in cos θ
  integer(I32) :: SysXc_Ylm_nTheta = 0

  !> Uniform points in φ (1: axial mode)
  integer(I32) :: SysXc_Ylm_nPhi = 0

  !> Gauss-Legendre nodes x_q = cos θ_q (SysXc_Ylm_nTheta; SysXc_Setup);
  !> the angles φ_k = 2π (k - 1)/nPhi
  real(R64), allocatable :: SysXc_Ylm_cosTheta(:)

  !> Axial mode: m blocks transformed separately, φ-averaged density
  logical :: SysXc_Ylm_axialQ = .false.

  !> Channel caps of the interaction layout (sysInteraction.ylm) of the
  !> Hartree source and potential
  integer(I32) :: SysXc_Ylm_lmaxSrc = 0
  integer(I32) :: SysXc_Ylm_mmaxSrc = 0

  !> Band limit of densities (2 lmax; |m| ≤ 2 mmax, 0 in the axial mode),
  !> the channels of the angular gradient
  integer(I32) :: SysXc_Ylm_lmaxDensity = 0
  integer(I32) :: SysXc_Ylm_mmaxDensity = 0

end module
