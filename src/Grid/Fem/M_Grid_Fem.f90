! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 3D Cartesian spectral-element (FEM) grid for molecules.
!>
!> Domain: the box [xmin, xmax] x [ymin, ymax] x [zmin, zmax]. The mesh is the
!> tensor product of three 1D element partitions, graded toward the "centers"
!> (nuclei): every center coordinate is an element boundary on every axis, so
!> every nucleus sits at an element vertex and the C0 basis captures its cusp.
!>
!> Basis: Lagrange polynomials on nLocals Gauss-Lobatto-Legendre (GLL) nodes
!> per element and axis (degree p = nLocals - 1), C0 continuous through the
!> shared element-boundary nodes. GLL quadrature makes the mass matrix
!> diagonal, W = Wx (x) Wy (x) Wz (lumped mass, as in FEDVR). Orbitals obey
!> homogeneous Dirichlet conditions on the box surface, so the two box-boundary
!> nodes of every axis are removed: n_a = nElements_a * (nLocals - 1) - 1
!> interior nodes per axis.
!>
!> Fields are nodal values f(ix, iy, iz), flattened with x fastest:
!>   i = ix + nx (iy - 1) + nx ny (iz - 1)          (Grid_Fem_Index)
!> and Grid_InnerProduct(f, g) = sum_i conjg(f_i) g_i w_i with the 3D product
!> weights w_i = Grid_Fem_weights(i).
!>
!> Operators. With the 1D stiffness K_a(i, j) = int L_i' L_j' dx (exact GLL
!> quadrature) the 3D stiffness is K3 = Kx(x)Wy(x)Wz + Wx(x)Ky(x)Wz + Wx(x)Wy(x)Kz
!> and the nodal Laplacian lap f = -W^-1 K3 f = -sum_a (W_a^-1 K_a) applied along
!> axis a; it is self-adjoint in the W metric. The fast diagonalization method
!> (FDM) diagonalizes every axis on its interior nodes,
!>   S_a^T K_a S_a = diag(lambda_a),   S_a^T W_a S_a = 1,
!> so with S3 = Sx (x) Sy (x) Sz and D(i, j, k) = lx(i) + ly(j) + lz(k):
!>   -lap = S3 D S3^T W3,   (K3 + s W3)^-1 = S3 (D + s)^-1 S3^T,
!> which gives exact shifted solves and any function of the Laplacian through
!> three dense mode products per axis (real DGEMMs).
!>
!> Electrostatics (Galerkin): -lap V = 4 pi rho is solved in weak form,
!>   K_II V_I = 4 pi load_I - K_IB V_B,   load_i = w_i rho_i,
!> with the box-surface values V_B from the free-space potential: a multipole
!> expansion of the load in real regular solid harmonics (Grid_Fem_SolvePoisson)
!> or the exact sum q/|r - R| (Grid_Fem_SolvePoissonPointCharges, Galerkin
!> point charges with load_i = q phi_i(R)).
!>
!> All module procedures are thread-safe (no module-level scratch) and may be
!> called inside OpenMP parallel regions; their own OpenMP loops are disabled
!> there (if (.not. omp_in_parallel())). The transforms, shifted solves and
!> Grid_Fem_SolvePoisson read their whole input before they write the output,
!> so they also work in place, although such a call violates the Fortran
!> argument aliasing rules (the code base does not rely on it).
!>
!> Code: three submodules, layered like their call graphs. S_Grid_Fem holds
!> the mesh (fabrication, Grid_Setup, inner product, node layout);
!> S_Grid_Fem_Laplacian the 1D operator tables of the axes (private), the
!> Laplacian, the FDM transforms and the shifted solves; its child
!> S_Grid_Fem_Poisson (a nested submodule) the Poisson solves, which reach
!> the operator tables and FDM kernels of the parent by host association.
module M_Grid_Fem
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! types
  !=============================================================================

  !> One Cartesian axis of the tensor-product mesh: element partition, interior
  !> GLL nodes and weights, and the spectrum of -d^2/dx^2 on them (the 1D
  !> operator tables of the 3D routines are private to S_Grid_Fem_Laplacian).
  !> Only interior nodes (box-boundary nodes removed, Dirichlet) carry degrees
  !> of freedom.
  type :: T_Grid_Fem_Axis
    !> Number of finite elements along the axis.
    integer(I32) :: nElements = 0
    !> Number of interior nodes, nElements * (nLocals - 1) - 1.
    integer(I32) :: nPoints = 0
    !> Box bounds of the axis (the two Dirichlet boundary nodes).
    real(R64) :: lowerBound = 0.0_R64
    real(R64) :: upperBound = 0.0_R64
    !> Element boundaries, bounds(0:nElements), bounds(0) = lowerBound,
    !> bounds(nElements) = upperBound.
    real(R64), allocatable :: bounds(:)
    !> Interior GLL nodes (nPoints); element boundaries are nodes, exactly.
    real(R64), allocatable :: points(:)
    !> Assembled GLL weights of the interior nodes (nPoints), the diagonal of W.
    real(R64), allocatable :: weights(:)
    !> Ascending eigenvalues of W^-1 K on the interior nodes (nPoints), i.e. of
    !> -d^2/dx^2 with Dirichlet conditions.
    real(R64), allocatable :: eigenvalues(:)
  end type

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure the FEM grid back-end: reads grid.fem, builds
    !> the element partitions of the three axes (so that Grid_nPoints is set on
    !> return) and binds Grid_InnerProduct and Grid_Setup. Grid_Setup builds
    !> the nodes, weights, FDM eigendecompositions (Grid_Fem_SetupAxes) and
    !> multipole tables (Grid_Fem_SetupPoisson).
    module subroutine Grid_Fem_Fabricate
    end subroutine

    !> Step of Grid_Setup (called by it only): nodes, weights and eigenvalues
    !> of the three axes and their 1D operator tables (stiffness band,
    !> boundary columns, FDM eigenvectors; private to S_Grid_Fem_Laplacian).
    !> Warns about sliver elements and stops if the smallest FDM eigenvalue of
    !> an axis is not positive.
    module subroutine Grid_Fem_SetupAxes
    end subroutine

    !> Step of Grid_Setup (called by it only, after Grid_Fem_SetupAxes and the
    !> mesh summary): tables of the multipole boundary values (solid-harmonic
    !> coefficients, face nodes, irregular solid harmonics there; private to
    !> S_Grid_Fem_Poisson) and Grid_Fem_multipoleRadius; prints the multipole
    !> lines of the mesh summary.
    module subroutine Grid_Fem_SetupPoisson
    end subroutine

    !> Nodal Laplacian with Dirichlet box conditions,
    !>   lapF = lap f = -sum_a (W_a^-1 K_a) f   (along axis a, banded stencil).
    module subroutine Grid_Fem_ApplyLaplacian(lapF, f)
      complex(R64), intent(out), contiguous :: lapF(:)
      complex(R64), intent(in), contiguous  :: f(:)
    end subroutine

    !> Coefficients in the eigenbasis of -lap: c = S3^T W3 f, so that
    !> f = S3 c and -lap f = S3 (D c) with D(i, j, k) = lx(i) + ly(j) + lz(k)
    !> (the flattened eigenvalue index follows the node layout, x fastest).
    module subroutine Grid_Fem_ToEigenbasis(c, f)
      complex(R64), intent(out), contiguous :: c(:)
      complex(R64), intent(in), contiguous  :: f(:)
    end subroutine

    !> Nodal values from eigenbasis coefficients: f = S3 c (inverse of
    !> Grid_Fem_ToEigenbasis).
    module subroutine Grid_Fem_FromEigenbasis(f, c)
      complex(R64), intent(out), contiguous :: f(:)
      complex(R64), intent(in), contiguous  :: c(:)
    end subroutine

    !> Galerkin solve of -lap pot = 4 pi rho for the weighted source
    !> load = w rho (as produced by an orbital pair product times the weights):
    !>   K_II pot = 4 pi load - K_IB V_B,
    !> with the box-surface values V_B from the multipole expansion of the load
    !> about Grid_Fem_multipoleCenter up to l = Grid_Fem_multipoleLmax. Real and
    !> imaginary parts are solved separately (the imaginary one only if nonzero);
    !> their moments come from one pass over the load. The expansion converges
    !> at the face nodes only for a load inside the sphere of radius
    !> Grid_Fem_multipoleRadius about the center: if more than 1e-4 of
    !> sum |load| lies outside, a warning is printed (once per run; the
    !> exponential tails of bound densities, ~1e-7 ... 1e-5, are harmless).
    module subroutine Grid_Fem_SolvePoisson(pot, load)
      complex(R64), intent(out), contiguous :: pot(:)
      complex(R64), intent(in), contiguous  :: load(:)
    end subroutine

    !> Galerkin potential pot ~ sum_I charges(I)/|r - positions(:, I)| of point
    !> charges strictly inside the box: load_i = sum_I q_I phi_i(R_I) (tensor
    !> Lagrange basis values; q_I delta_iJ if R_I is the node J), exact
    !> box-surface values sum_I q_I/|r_B - R_I|. The value at a charge sitting on
    !> a node is finite (roughly 30 ... 50 q/hMin for nLocals 5-6, depending on
    !> the mesh, about 85 q/hMin for nLocals 8; hMin the local
    !> grid.fem.centers.hMin).
    module subroutine Grid_Fem_SolvePoissonPointCharges(pot, charges, positions)
      complex(R64), intent(out), contiguous :: pot(:)
      real(R64), intent(in), contiguous     :: charges(:)
      real(R64), intent(in), contiguous     :: positions(:, :)
    end subroutine

    !> Flattened index of the interior node (ix, iy, iz):
    !>   i = ix + nx (iy - 1) + nx ny (iz - 1).
    pure module function Grid_Fem_Index(ix, iy, iz) result(i)
      integer(I32), intent(in) :: ix
      integer(I32), intent(in) :: iy
      integer(I32), intent(in) :: iz
      integer(I32) :: i
    end function

    !> Node coordinate along a direction, coord(i) = direction . r_i at every
    !> interior node i (allocated to Grid_nPoints; r_i measured from the
    !> coordinate origin; after Grid_Setup), evaluated as
    !> direction(1) x + (direction(2) y + direction(3) z): a unit axis vector
    !> gives the node coordinates x, y or z exactly, a unit vector e the
    !> coordinate along e (e.g. e_pol . r of the length gauge, kicks, dipoles).
    module subroutine Grid_Fem_FillNodeCoordinate(coord, direction)
      real(R64), allocatable, intent(out) :: coord(:)
      real(R64), intent(in)               :: direction(3)
    end subroutine
  end interface

  !> Exact shifted solve u = (-lap + shift)^-1 f = S3 (D + shift)^-1 S3^T W3 f
  !> for a real or a complex shift (generic). The mode products stay real;
  !> only the diagonal (D + shift)^-1 is complex for a complex shift, where
  !> the real and imaginary parts of the coefficients mix.
  interface Grid_Fem_SolveShiftedLaplacian
    !> Real shift; requires shift > -(smallest eigenvalue of -lap). Real and
    !> imaginary parts of f are solved separately (the imaginary one only if
    !> nonzero).
    module subroutine Grid_Fem_SolveShiftedLaplacianReal(u, f, shift)
      complex(R64), intent(out), contiguous :: u(:)
      complex(R64), intent(in), contiguous  :: f(:)
      real(R64), intent(in)                 :: shift
    end subroutine

    !> Complex shift, e.g. 2/(i dt) for the Crank-Nicolson factor
    !> 1 + i dt (-lap)/2 = (i dt/2) (-lap + 2/(i dt)). A shift with zero
    !> imaginary part takes the real-shift path (and its validity check);
    !> otherwise D + shift never vanishes (D is real).
    module subroutine Grid_Fem_SolveShiftedLaplacianComplex(u, f, shift)
      complex(R64), intent(out), contiguous :: u(:)
      complex(R64), intent(in), contiguous  :: f(:)
      complex(R64), intent(in)              :: shift
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Number of GLL nodes per element edge (polynomial degree nLocals - 1).
  integer(I32) :: Grid_Fem_nLocals = 0

  !> The x, y and z axes of the tensor-product mesh.
  type(T_Grid_Fem_Axis) :: Grid_Fem_axes(3)

  !> Interior nodes per axis (copies of Grid_Fem_axes(:) % nPoints, set at
  !> fabrication); Grid_nPoints = nx ny nz.
  integer(I32) :: Grid_Fem_nx = 0
  integer(I32) :: Grid_Fem_ny = 0
  integer(I32) :: Grid_Fem_nz = 0

  !> Number of centers (nuclei) the mesh is graded toward.
  integer(I32) :: Grid_Fem_nCenters = 0

  !> Center positions (3, nCenters); every coordinate is an element boundary
  !> of the corresponding axis unless explicit bounds override the grading or
  !> it was merged with a closer than 0.25 min(hMin) coordinate of another
  !> center (a note is printed at fabrication).
  real(R64), allocatable :: Grid_Fem_centers(:, :)

  !> Smallest element size at each center (nCenters).
  real(R64), allocatable :: Grid_Fem_centerHMin(:)

  !> 3D product quadrature weights w = wx(ix) wy(iy) wz(iz) (Grid_nPoints).
  real(R64), allocatable :: Grid_Fem_weights(:)

  !> Largest l of the multipole expansion of the Poisson boundary values.
  integer(I32) :: Grid_Fem_multipoleLmax = 0

  !> Expansion center of the Poisson boundary values (default: the centroid of
  !> the centers, the box center without centers).
  real(R64) :: Grid_Fem_multipoleCenter(3) = 0.0_R64

  !> Convergence radius of the multipole boundary values: the distance from
  !> Grid_Fem_multipoleCenter to the nearest box-surface node (set in
  !> Grid_Setup). The expansion converges there only for a load inside this
  !> sphere; elongated boxes or charge far from the center violate it.
  real(R64) :: Grid_Fem_multipoleRadius = 0.0_R64

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
