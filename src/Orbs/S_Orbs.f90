! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Implementation submodule for M_Orbs.
!>
!> This submodule provides the concrete implementations for orbital operations:
!> - Fabrication: procedure pointer binding and configuration parsing
!> - Orthonormalization: per-body-type Gram-Schmidt via Grid backend
!> - Subspace projection: gauge enforcement via Grid backend
!> - Persistence: binary file I/O for orbital coefficients
!>
!> Design rationale
!> ----------------
!> Orbital operations are thin wrappers around Grid operations. The key insight
!> is that each body type's orbitals form an independent subspace and must be
!> orthonormalized separately. This ensures that fermions/bosons of different
!> species remain distinguishable and that the many-body wavefunction retains
!> proper antisymmetry/symmetry structure.
!>
!> The loop pattern `do ibt = 1, nBodyTypes` appears in all operations and
!> partitions the orbital matrix into body-type slices using the indexing:
!>   startOrb = nOrbsStart(ibt), endOrb = nOrbsEnd(ibt)
!>
!> Dependencies
!> ------------
!> - M_Grid: provides metric-aware orthonormalization and projection
!> - M_Method_Mb: defines nBodyTypes (number of particle species)
!> - M_Method_Mb_OrbBased: defines nOrbs(ibt), nOrbsSum (orbital partitioning)
!> - M_Utils_Json: JSON configuration access
!> - M_Utils_DataStorage: binary I/O for SaveOrbs
submodule(M_Orbs) S_Orbs

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Bind procedure pointers and configure the orbital representation.
!>
!> Configuration keys (JSON)
!> -------------------------
!> - `orbs.restrictedQ` (logical, default: false):
!>     If true, use spin-restricted representation where a single set of
!>     spatial orbitals is shared between spin-up and spin-down electrons.
!>     This halves `Orbs_nOrbsInState` and enables Hamiltonian symmetries.
!>
!> Procedure pointer bindings
!> --------------------------
!> - Orbs_ProjectOnSubspace => ProjectOnSubspace
!> - Orbs_Orthonormalize    => Orthonormalize
!> - Orbs_SaveOrbs          => SaveOrbs
!>
!> Post-conditions
!> ---------------
!> - `Orbs_restrictedQ` is set from JSON configuration
!> - `Orbs_nOrbsInState` is computed:
!>     * Restricted:   nOrbsSum / 2 (shared spatial orbitals)
!>     * Unrestricted: nOrbsSum     (independent orbitals per body type)
  module subroutine Orbs_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    call Say_Fabricate("orbs")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Orbs_ProjectOnSubspace => ProjectOnSubspace
    Orbs_Orthonormalize => Orthonormalize
    Orbs_SaveOrbs => SaveOrbs

    Orbs_restrictedQ = Json_Get("orbs.restrictedQ", .false.)

    ! The restricted representation stores a single spatial orbital set shared
    ! by spin up and spin down; it is only meaningful for exactly two body
    ! types (the spins) with identical orbital counts
    if (Orbs_restrictedQ) then
      if (Method_Mb_nBodyTypes .ne. 2) then
        error stop "orbs.restrictedQ requires exactly 2 body types (spin up/down)"
      end if
      if (Method_Mb_OrbBased_nOrbs(1) .ne. Method_Mb_OrbBased_nOrbs(2)) then
        error stop "orbs.restrictedQ requires equal nOrbs for both body types"
      end if
    end if

    if (Orbs_restrictedQ) Orbs_nOrbsInState = Method_Mb_OrbBased_nOrbsSum / 2
    if (.not. Orbs_restrictedQ) Orbs_nOrbsInState = Method_Mb_OrbBased_nOrbsSum

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Gauge-align orbitals with a reference set.
!>
!> The overlap matrix M = <orbs_i | refOrbs_j> is decomposed as
!> M = U * S * V^H, then the unitary W = U * V^H is applied:
!>   orbs <- orbs * W
!> This gives the orthonormal basis in span(orbs) closest to refOrbs.
  module subroutine Orbs_AlignOnReference(orbs, refOrbs)
    use M_Grid
    use M_Utils_LapackLib

    complex(R64), intent(inout), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: refOrbs(:, :)

    complex(R64), allocatable :: overlap(:, :), u(:, :), vt(:, :)
    real(R64), allocatable :: singularVals(:)
    integer(I32) :: nOrbs, iOrb, jOrb

    nOrbs = size(orbs, 2)

    allocate (overlap(nOrbs, nOrbs))
    do jOrb = 1, nOrbs
      do iOrb = 1, nOrbs
        overlap(iOrb, jOrb) = Grid_InnerProduct(orbs(:, iOrb), refOrbs(:, jOrb))
      end do
    end do

    call LapackLib_Svd(u, singularVals, vt, overlap)
    orbs = matmul(orbs, matmul(u, vt))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Orthonormalize orbital columns independently for each body type.
!>
!> Algorithm
!> ---------
!> For each body type ibt in [1, nBodyTypes]:
!>   1. Extract the orbital slice orbs(:, startOrb:endOrb) for this body type
!>   2. Apply Grid_Orthonormalize which performs modified Gram-Schmidt with
!>      the grid's metric tensor (handles non-uniform grids, FEDVR, etc.)
!>   3. Result: <φ_i^(bt) | φ_j^(bt)> = δ_ij within body type
!>
!> Why per-body-type?
!> ------------------
!> Cross-body-type orthogonality is not required (and would be unphysical for
!> distinguishable particles). Each species maintains its own orthonormal
!> orbital basis independently.
!>
!> Complexity: O(nBodyTypes * nOrbs(ibt)^2 * nBasis) dominated by inner products.
  subroutine Orthonormalize(orbs)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Grid

    complex(R64), intent(inout), contiguous  :: orbs(:, :)

    integer(I32) :: ibt, startOrb, endOrb

    ! Restricted: the stored columns are the single spatial set shared by both
    ! spins; orthonormalize them as one block (the body-type partition refers
    ! to the full spin-orbital set, which is not stored)
    if (Orbs_restrictedQ) then
      call Grid_Orthonormalize(orbs)
      return
    end if

    startOrb = 1
    do ibt = 1, Method_Mb_nBodyTypes

      endOrb = startOrb + Method_Mb_OrbBased_nOrbs(ibt) - 1

      call Grid_Orthonormalize(orbs(:, startOrb:endOrb))

      startOrb = endOrb + 1

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Project orbital derivatives onto the orthogonal complement of the orbital
!> subspace, independently for each body type.
!>
!> Mathematical operation
!> ----------------------
!> For each column ψ in dOrbs and the reference orbitals {φ_j} in orbs:
!>   ψ' = ψ - Σ_j <φ_j|ψ> φ_j
!>
!> This removes the component of ψ parallel to span{φ_1, ..., φ_n}.
!>
!> Physical motivation (MCTDH gauge)
!> ---------------------------------
!> In MCTDH/MCTDHF/MCTDHB propagation, the orbital time derivatives dφ/dt
!> must be orthogonal to the current orbital space to avoid redundant
!> parametrization (gauge freedom). This projection enforces the standard
!> "parallel transport" gauge condition: <φ_j | dφ_i/dt> = 0.
!>
!> Why per-body-type?
!> ------------------
!> The gauge condition applies within each body type independently. Orbitals
!> of different species are already distinguishable and need not (and should
!> not) be projected against each other.
  subroutine ProjectOnSubspace(dOrbs, orbs)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Grid

    complex(R64), intent(inout), contiguous  :: dOrbs(:, :)
    complex(R64), intent(in), contiguous     :: orbs(:, :)

    integer(I32) :: ibt, startOrb, endOrb

    ! Restricted: single shared spatial set, project as one block
    if (Orbs_restrictedQ) then
      call Grid_ProjectOnSubspace(dOrbs, orbs)
      return
    end if

    startOrb = 1
    do ibt = 1, Method_Mb_nBodyTypes

      endOrb = startOrb + Method_Mb_OrbBased_nOrbs(ibt) - 1

      call Grid_ProjectOnSubspace(dOrbs(:, startOrb:endOrb), orbs(:, startOrb:endOrb))

      startOrb = endOrb + 1

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Persist orbital coefficients to binary files, one file per orbital.
!>
!> File naming convention
!> ----------------------
!> Each orbital is saved to: `orbBB_II.in`
!>   - BB: two-digit body type index (01, 02, ...)
!>   - II: two-digit orbital index within that body type (01, 02, ...)
!>
!> Example for 2 body types with 3 and 2 orbitals respectively:
!>   orb01_01.in, orb01_02.in, orb01_03.in  (body type 1)
!>   orb02_01.in, orb02_02.in               (body type 2)
!>
!> File format
!> -----------
!> Raw binary via M_Utils_DataStorage. Each file contains the complex(R64)
!> coefficient vector of length Grid_nPoints. Files can be reloaded via
!> OrbsInit with backend "Load".
!>
!> Use cases
!> ---------
!> - Checkpoint/restart for long propagations
!> - Post-processing and visualization
!> - Initial guess for subsequent calculations
  subroutine SaveOrbs(orbs)
    use M_Utils_DataStorage
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(in), contiguous  :: orbs(:, :)

    integer(I32) :: index, ibt, i1, nBt
    character(len=256) :: filename

    ! Restricted: only the shared spatial set exists, save it as body type 1
    nBt = Method_Mb_nBodyTypes
    if (Orbs_restrictedQ) nBt = 1

    i1 = 0
    do ibt = 1, nBt
      do index = 1, Method_Mb_OrbBased_nOrbs(ibt)
        i1 = i1 + 1

        write (filename, '("orb",I2.2,"_",I2.2,".in")') ibt, index
        call SaveData(orbs(:, i1), trim(filename), storage_size(orbs(:, i1)))

      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Mix occupied orbitals using density matrix formulation.
  module subroutine Orbs_MixOccupiedSpace(orbsOld, evecs, alpha, orbsOut, lambdaDiscarded)
    use M_Grid
    use M_Utils_LapackLib

    complex(R64), intent(in), contiguous  :: orbsOld(:,:), evecs(:,:)
    real(R64),    intent(in)  :: alpha
    complex(R64), intent(out), contiguous :: orbsOut(:,:)
    real(R64),    intent(out), allocatable :: lambdaDiscarded(:)

    integer(I32) :: n, i, j, k
    complex(R64), allocatable :: A(:,:), G(:,:), G_evecs(:,:)
    real(R64),    allocatable :: lambda(:)

    n = size(orbsOld, 2)
    allocate (A(size(orbsOld, 1), 2 * n))
    allocate (G(2 * n, 2 * n))

    A(:, 1:n) = sqrt(1.0_R64 - alpha) * orbsOld
    A(:, n + 1:2 * n) = sqrt(alpha) * evecs

    do j = 1, 2 * n
      do i = 1, 2 * n
        G(i, j) = Grid_InnerProduct(A(:, i), A(:, j))
      end do
    end do

    call LapackLib_DiagonalizeGeneric(lambda, G_evecs, G, .true.)

    orbsOut = 0.0_R64
    do i = 1, n
      k = n + i
      do j = 1, 2 * n
        orbsOut(:, i) = orbsOut(:, i) + G_evecs(j, k) * A(:, j)
      end do
      orbsOut(:, i) = orbsOut(:, i) / sqrt(lambda(k))
    end do

    lambdaDiscarded = lambda(1:n)

  end subroutine

end submodule
