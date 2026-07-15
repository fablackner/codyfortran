! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Prolate-spheroidal grid API for homonuclear diatomic molecules.
!>
!> Two-center coordinate system with the nuclei at the foci z = -a and z = +a,
!> where a = R/2 and R is the (fixed) internuclear distance:
!>
!>   xi  in [1, ximax]   ("radial-like"),  r1 = a(xi+eta), r2 = a(xi-eta)
!>   eta in [-1, 1]      ("angular-like"), z = a*xi*eta
!>   phi in [0, 2pi)     (azimuthal angle around the molecular axis)
!>
!> The volume element is dV = a^3 (xi^2 - eta^2) dxi deta dphi, which cancels
!> both nuclear Coulomb singularities analytically. Fields are stored in
!> azimuthal channels (m is a good quantum number of the linear molecule):
!>
!>   psi(xi, eta, phi) = sum_m f_m(xi, eta) e^(i m phi) / sqrt(2 pi)
!>
!> Flattened layout (mirrors the Ylm channel-outer convention): the m channel
!> runs outermost in the interleaved order m = 0, -1, +1, -2, +2, ..., then
!> eta, then xi fastest, so each m channel is one contiguous block of
!> nXi*nEta spatial values. The interleaved order makes a channel's block
!> position independent of the field's own mmax, so fields with different
!> mmax (orbitals vs. interaction potentials) share the same accessors.
!>
!> Discretization: FEDVR in xi (back-end selected at fabrication; xi = 1 is a
!> grid point since sigma orbitals are nonzero on the internuclear axis, ximax
!> carries a Dirichlet condition) and a single-element Gauss-Legendre DVR in
!> eta (interior nodes only, so the endpoint factors (1 - eta^2) never vanish
!> on the grid).
!>
!> @note For odd |m| the physical fields carry a factor ((xi^2-1)(1-eta^2))^(|m|/2)
!> that is not polynomial in (xi, eta); the polynomial DVR basis then converges
!> only algebraically. Even-m channels (sigma, delta, ...) are fully spectral.
module M_Grid_Prolate
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a prolate-spheroidal grid back-end.
    module subroutine Grid_Prolate_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Internuclear distance R (the foci sit at z = -R/2 and z = +R/2).
  real(R64) :: Grid_Prolate_R

  !> Half the internuclear distance, a = R/2 (focal distance).
  real(R64) :: Grid_Prolate_a

  !> Outer boundary of the xi coordinate (box radius along z is a*ximax).
  real(R64) :: Grid_Prolate_ximax

  !> Maximum azimuthal quantum number; channels run m = -mmax..mmax.
  integer(I32) :: Grid_Prolate_mmax

  !> Number of xi grid points (filled at fabrication by the xi back-end).
  integer(I32) :: Grid_Prolate_nXi

  !> Number of eta grid points (Gauss-Legendre, interior only).
  integer(I32) :: Grid_Prolate_nEta

  !> Number of spatial points per m channel, nXi * nEta.
  integer(I32) :: Grid_Prolate_nSpatial

  !> xi quadrature nodes (xi = 1 included, ximax excluded).
  real(R64), allocatable :: Grid_Prolate_xiPoints(:)

  !> xi quadrature weights (bare, without the (xi^2 - eta^2) Jacobian).
  real(R64), allocatable :: Grid_Prolate_xiWeights(:)

  !> eta quadrature nodes (Gauss-Legendre on [-1, 1], endpoints excluded).
  real(R64), allocatable :: Grid_Prolate_etaPoints(:)

  !> eta quadrature weights.
  real(R64), allocatable :: Grid_Prolate_etaWeights(:)

  !> Symmetric Sturm-Liouville DVR matrix of the xi kinetic term,
  !>   S_xi(i,j) = sum_k w_k (xi_k^2 - 1) D(k,i) D(k,j)
  !> so that [d/dxi (xi^2-1) d/dxi f]_i = -(1/w_i) sum_j S_xi(i,j) f_j.
  !> Filled by the xi back-end during setup.
  real(R64), allocatable :: Grid_Prolate_xiKinMatrix(:, :)

  !> Symmetric Sturm-Liouville DVR matrix of the eta kinetic term,
  !>   S_eta(i,j) = sum_k w_k (1 - eta_k^2) D(k,i) D(k,j)
  !> so that [d/deta (1-eta^2) d/deta f]_j = -(1/w_j) sum_j' S_eta(j,j') f_j'.
  real(R64), allocatable :: Grid_Prolate_etaKinMatrix(:, :)

  !> Spatial metric weights per (xi, eta) point in channel-block order
  !> (xi fastest): w_xi * w_eta * a^3 * (xi^2 - eta^2). One channel's worth;
  !> identical for every m block.
  real(R64), allocatable :: Grid_Prolate_spatialWeights(:)

  !> Total integration weights for each (xi, eta, m) entry in flattened order
  !> (spatialWeights repeated per m block).
  real(R64), allocatable :: Grid_Prolate_weights(:)

  !> Component-wise coordinates for each flattened entry.
  real(R64), allocatable :: Grid_Prolate_xiCoord(:)
  real(R64), allocatable :: Grid_Prolate_etaCoord(:)
  integer(I32), allocatable :: Grid_Prolate_mCoord(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the spatial-product routine producing azimuthal channels.
  procedure(I_Grid_Prolate_SpatialProduct), pointer :: Grid_Prolate_SpatialProduct
  abstract interface
    !> Compute fOut = fConjg^(*) * f expanded in e^(i m phi)/sqrt(2 pi) channels.
    !>
    !> The (xi, eta) dependence multiplies pointwise; the azimuthal channels
    !> convolve, picking up one factor 1/sqrt(2 pi) from the channel basis:
    !>   fOut_m = (1/sqrt(2 pi)) sum_(m1) fConjg_(m1) f_(m1+m)      (conjgQ_)
    !>   fOut_m = (1/sqrt(2 pi)) sum_(m1) fConjg_(m1) f_(m-m1)      (otherwise)
    !> Optional flags:
    !>  - conjgQ_: when true, conjugates fConjg (bra factor).
    !>  - withWeightsQ_: when true, multiplies the result by the spatial
    !>    metric weights (quadrature x Jacobian).
    pure subroutine I_Grid_Prolate_SpatialProduct(fOut, fConjg, f, mMaxOut, mMaxConjg, mMax, conjgQ_, withWeightsQ_)
      import :: I32, R64
      complex(R64), intent(out), contiguous :: fOut(:)
      complex(R64), intent(in), contiguous  :: fConjg(:)
      complex(R64), intent(in), contiguous  :: f(:)
      integer(I32), intent(in)              :: mMaxOut
      integer(I32), intent(in)              :: mMaxConjg
      integer(I32), intent(in)              :: mMax
      logical, intent(in), optional         :: conjgQ_
      logical, intent(in), optional         :: withWeightsQ_
    end subroutine
  end interface

  !> Pointer to the setup procedure of the xi back-end (nodes, weights,
  !> Sturm-Liouville matrix).
  procedure(I_Grid_Prolate_Setup), pointer :: Grid_Prolate_Setup
  abstract interface
    !> Initialize the active xi back-end.
    subroutine I_Grid_Prolate_Setup
    end subroutine
  end interface

  !> Pointer to extract a single m channel from a function on the grid.
  procedure(I_Grid_Prolate_GetMComponent), pointer :: Grid_Prolate_GetMComponent
  abstract interface
    !> Extract the spatial function f_m(xi, eta) from a full field. The output
    !> is the contiguous channel block of size nSpatial (xi fastest).
    pure subroutine I_Grid_Prolate_GetMComponent(fMOut, m, f)
      import :: I32, R64
      !> Output spatial function for the requested m channel.
      complex(R64), intent(out) :: fMOut(:)
      !> The azimuthal quantum number (-mmax <= m <= mmax of the field f).
      integer(I32), intent(in) :: m
      !> The full function on the prolate grid (flattened order).
      complex(R64), intent(in) :: f(:)
    end subroutine
  end interface

  !> Pointer to set/overwrite a single m channel in a function on the grid.
  procedure(I_Grid_Prolate_SetMComponent), pointer :: Grid_Prolate_SetMComponent
  abstract interface
    !> Set the m channel of fOut to the provided spatial function.
    pure subroutine I_Grid_Prolate_SetMComponent(fOut, m, fM)
      import :: I32, R64
      !> The full function on the prolate grid after setting the channel.
      complex(R64), intent(out) :: fOut(:)
      !> The azimuthal quantum number.
      integer(I32), intent(in) :: m
      !> The input spatial function for the m channel.
      complex(R64), intent(in) :: fM(:)
    end subroutine
  end interface

  !> Pointer to add/accumulate a single m channel into a function on the grid.
  procedure(I_Grid_Prolate_AddMComponent), pointer :: Grid_Prolate_AddMComponent
  abstract interface
    !> Add the provided spatial function into the m channel of fInOut.
    pure subroutine I_Grid_Prolate_AddMComponent(fInOut, m, fM)
      import :: I32, R64
      !> The full function on the prolate grid to which the channel is added.
      complex(R64), intent(inout) :: fInOut(:)
      !> The azimuthal quantum number.
      integer(I32), intent(in) :: m
      !> The input spatial function for the m channel.
      complex(R64), intent(in) :: fM(:)
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
