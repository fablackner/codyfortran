! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Two-body RDM tensor algebra: symmetrization, antisymmetrization, and spin-block decomposition.
!>
!> @details
!> This module provides fundamental operations for manipulating 4-index tensors
!> (rank-4) that arise in 2-RDM theory. The key operations include:
!>
!>   - **Symmetrization/Antisymmetrization**: Project tensors onto symmetric or
!>     antisymmetric subspaces under particle exchange (fermionic/bosonic).
!>
!>   - **Wedge Products**: Construct 2-body tensors as (anti)symmetrized products
!>     of two 1-body density matrices (outer products with exchange terms).
!>
!>   - **Spin-Block Decomposition**: Decompose the α-β 2-RDM into singlet (S=0)
!>     and triplet (S=1) irreducible representations under spin SU(2).
!>
!> ## Tensor Indexing Convention
!>
!> All routines follow the physicist's convention for RDM indices:
!>
!>     A2(i1, i2, j1, j2) ≡ ⟨Ψ| a†_{i1} a†_{i2} a_{j2} a_{j1} |Ψ⟩
!>                        ≡ ⟨i1, i2 | Â | j1, j2⟩
!>
!> where:
!>   - (i1, i2) are **creation** (row) indices, mapping to bra orbitals
!>   - (j1, j2) are **annihilation** (column) indices, mapping to ket orbitals
!>   - Index 1 refers to particle 1, index 2 to particle 2
!>
!> ## Symmetry Properties
!>
!> For **fermionic** systems, the 2-RDM satisfies antisymmetry under exchange:
!>
!>     D²(i1,i2,j1,j2) = -D²(i2,i1,j1,j2)    [exchange of creation indices]
!>     D²(i1,i2,j1,j2) = -D²(i1,i2,j2,j1)    [exchange of annihilation indices]
!>     D²(i1,i2,j1,j2) = +D²(j1,j2,i1,i2)*   [Hermitian conjugate]
!>
!> For **bosonic** systems, the 2-RDM is fully symmetric under particle exchange.
!>
!> ## Related Modules
!>
!>   - `M_Utils_RdmCumulants`: Cumulant decomposition D² = D¹∧D¹ + Δ²
!>   - `M_Utils_RdmHole`: Q-RDM (hole density) and G-RDM construction
!>   - `M_Utils_RdmTrace`: Partial trace operations
!>   - `M_Utils_RdmPacking`: Pack/unpack to state vectors
!>
!> ## References
!>
!>   - Coleman, A.J. "Structure of Fermion Density Matrices", Rev. Mod. Phys. 35, 668 (1963)
!>   - Kutzelnigg & Mukherjee, J. Chem. Phys. 110, 2800 (1999)
!>
!> @note All routines assume complex(R64) data type for compatibility with
!>       time-dependent simulations where RDMs may become complex.
!>
!> @see AGENTS.md for developer onboarding and architectural context.
module M_Utils_RdmAlgebra
  use M_Utils_Types

  implicit none

  private
  public :: RdmAlgebra_Antisymmetrize
  public :: RdmAlgebra_Symmetrize
  public :: RdmAlgebra_FillAntisymProduct
  public :: RdmAlgebra_FillSymProduct
  public :: RdmAlgebra_FillSingletPart
  public :: RdmAlgebra_FillTripletPart

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Antisymmetrize a rank-4 tensor under particle exchange (fermionic projection).
!>
!> @details
!> Projects A2 onto the antisymmetric subspace satisfying fermionic exchange:
!>
!>     A2(i1,i2,j1,j2) = -A2(i2,i1,j1,j2) = -A2(i1,i2,j2,j1) = +A2(i2,i1,j2,j1)
!>
!> The projection operator P_A (idempotent: P_A² = P_A) is applied:
!>
!>     A2 ← P_A[A2] = (1/4) * [A2(i1,i2,j1,j2) - A2(i2,i1,j1,j2)
!>                                             - A2(i1,i2,j2,j1) + A2(i2,i1,j2,j1)]
!>
!> **Use Cases:**
!>   - Enforcing fermionic antisymmetry on approximate 2-RDMs
!>   - Post-processing numerical noise that breaks exchange symmetry
!>   - Validating that a tensor lies in the antisymmetric subspace
!>
!> **Complexity:** O(nO⁴) with temporary allocation of size nO⁴.
!>
!> @param[inout] A2  Rank-4 tensor of shape (nO, nO, nO, nO), modified in-place.
!>                   After call, A2 satisfies fermionic exchange antisymmetry.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmAlgebra_Antisymmetrize(A2)

    complex(R64), intent(inout), contiguous :: A2(:, :, :, :)

    integer(I32) :: nO
    integer(I32) :: i1, i2, j1, j2
    complex(R64), allocatable :: A2Tmp(:, :, :, :)

    nO = size(A2, 1)
    allocate (A2Tmp, source=A2)

    do j2 = 1, nO
      do j1 = 1, nO

        do i2 = 1, nO
          do i1 = 1, nO

            A2(i1, i2, j1, j2) = 0.25_R64 * (A2Tmp(i1, i2, j1, j2) - &
                                             A2Tmp(i2, i1, j1, j2) - &
                                             A2Tmp(i1, i2, j2, j1) + &
                                             A2Tmp(i2, i1, j2, j1))

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Symmetrize a rank-4 tensor under particle exchange (bosonic projection).
!>
!> @details
!> Projects A2 onto the symmetric subspace satisfying bosonic exchange:
!>
!>     A2(i1,i2,j1,j2) = +A2(i2,i1,j1,j2) = +A2(i1,i2,j2,j1) = +A2(i2,i1,j2,j1)
!>
!> The projection operator P_S (idempotent: P_S² = P_S) is applied:
!>
!>     A2 ← P_S[A2] = (1/4) * [A2(i1,i2,j1,j2) + A2(i2,i1,j1,j2)
!>                                             + A2(i1,i2,j2,j1) + A2(i2,i1,j2,j1)]
!>
!> **Use Cases:**
!>   - Enforcing bosonic symmetry on approximate 2-RDMs
!>   - Constructing symmetric components in spin-block decomposition
!>   - Post-processing numerical noise that breaks exchange symmetry
!>
!> **Complexity:** O(nO⁴) with temporary allocation of size nO⁴.
!>
!> @param[inout] A2  Rank-4 tensor of shape (nO, nO, nO, nO), modified in-place.
!>                   After call, A2 satisfies bosonic exchange symmetry.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmAlgebra_Symmetrize(A2)

    complex(R64), intent(inout), contiguous :: A2(:, :, :, :)

    integer(I32) :: nO
    integer(I32) :: i1, i2, j1, j2
    complex(R64), allocatable :: A2Tmp(:, :, :, :)

    nO = size(A2, 1)
    allocate (A2Tmp, source=A2)

    do j2 = 1, nO
      do j1 = 1, nO

        do i2 = 1, nO
          do i1 = 1, nO

            A2(i1, i2, j1, j2) = 0.25_R64 * (A2Tmp(i1, i2, j1, j2) + &
                                             A2Tmp(i2, i1, j1, j2) + &
                                             A2Tmp(i1, i2, j2, j1) + &
                                             A2Tmp(i2, i1, j2, j1))

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Construct antisymmetric wedge product of two 1-body matrices.
!>
!> @details
!> Computes the exterior (wedge) product A0 ∧ A1 as the antisymmetrized
!> outer product of two 1-body matrices:
!>
!>     A2(i1,i2,j1,j2) = (1/4) * [ A0(i1,j1)*A1(i2,j2) - A0(i2,j1)*A1(i1,j2)
!>                               - A0(i1,j2)*A1(i2,j1) + A0(i2,j2)*A1(i1,j1) ]
!>
!> This corresponds to the **fermionic** product with antisymmetrization,
!> implementing the Grassmann algebra wedge product for creation operators.
!>
!> **Physical Interpretation:**
!> When A0 = A1 = D¹ (the 1-RDM), the result is the "free-particle" or
!> "disconnected" 2-RDM representing uncorrelated two-body density:
!>
!>     D²_free = D¹ ∧ D¹
!>
!> This is the mean-field (Hartree-Fock) approximation. The difference
!> D² - D²_free defines the 2-body cumulant Δ² measuring correlation.
!>
!> **Identity:** For self-product: A ∧ A = 0 when A has rank 1.
!>
!> **Complexity:** O(nO⁴) with O(nO⁴) output allocation.
!>
!> @param[out]    A2  Rank-4 antisymmetric tensor (allocated if needed), shape (nO, nO, nO, nO).
!> @param[in]     A0  First 1-body matrix, shape (nO, nO).
!> @param[in]     A1  Second 1-body matrix, shape (nO, nO).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmAlgebra_FillAntisymProduct(A2, A0, A1)

    complex(R64), intent(out), allocatable :: A2(:, :, :, :)
    complex(R64), intent(in), contiguous                :: A0(:, :)
    complex(R64), intent(in), contiguous                :: A1(:, :)

    integer(I32) :: nO
    integer(I32) :: i1, i2, j1, j2

    nO = size(A0, 1)

    if (.not. allocated(A2)) allocate (A2(nO, nO, nO, nO))

    Do i1 = 1, nO
      Do j1 = 1, nO
        Do i2 = 1, nO
          Do j2 = 1, nO

            A2(i1, i2, j1, j2) = 0.25_R64 * (A0(i1, j1) * A1(i2, j2) - &
                                             A0(i2, j1) * A1(i1, j2) - &
                                             A0(i1, j2) * A1(i2, j1) + &
                                             A0(i2, j2) * A1(i1, j1))

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Construct symmetric product of two 1-body matrices.
!>
!> @details
!> Computes the symmetrized outer product A0 ⊙ A1 of two 1-body matrices:
!>
!>     A2(i1,i2,j1,j2) = (1/4) * [ A0(i1,j1)*A1(i2,j2) + A0(i2,j1)*A1(i1,j2)
!>                               + A0(i1,j2)*A1(i2,j1) + A0(i2,j2)*A1(i1,j1) ]
!>
!> This corresponds to the **bosonic** product with symmetrization,
!> analogous to the symmetric tensor product for bosonic creation operators.
!>
!> **Physical Interpretation:**
!> For bosonic systems, this gives the uncorrelated 2-body density.
!> For spin-singlet fermions, this structure appears in the singlet
!> component of the α-β 2-RDM decomposition.
!>
!> **Complexity:** O(nO⁴) with O(nO⁴) output allocation.
!>
!> @param[out]    A2  Rank-4 symmetric tensor (allocated if needed), shape (nO, nO, nO, nO).
!> @param[in]     A0  First 1-body matrix, shape (nO, nO).
!> @param[in]     A1  Second 1-body matrix, shape (nO, nO).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmAlgebra_FillSymProduct(A2, A0, A1)

    complex(R64), intent(out), allocatable :: A2(:, :, :, :)
    complex(R64), intent(in), contiguous                :: A0(:, :)
    complex(R64), intent(in), contiguous                :: A1(:, :)

    integer(I32) :: nO
    integer(I32) :: i1, i2, j1, j2

    nO = size(A0, 1)

    if (.not. allocated(A2)) allocate (A2(nO, nO, nO, nO))

    Do i1 = 1, nO
      Do j1 = 1, nO
        Do i2 = 1, nO
          Do j2 = 1, nO

            A2(i1, i2, j1, j2) = 0.25_R64 * (A0(i1, j1) * A1(i2, j2) + &
                                             A0(i2, j1) * A1(i1, j2) + &
                                             A0(i1, j2) * A1(i2, j1) + &
                                             A0(i2, j2) * A1(i1, j1))

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Extract singlet (S=0) component from α-β 2-RDM.
!>
!> @details
!> Decomposes the spin-resolved α-β 2-RDM into spin irreducible representations
!> (irreps) under the SU(2) spin rotation group. The singlet (S=0) component
!> is **symmetric** under combined exchange of spatial and spin indices:
!>
!>     D²_singlet(i1,i2,j1,j2) = (1/4) * [ D²_αβ(i1,i2,j1,j2) + D²_αβ(i2,i1,j2,j1)
!>                                       + D²_αβ(i1,i2,j2,j1) + D²_αβ(i2,i1,j1,j2) ]
!>
!> ## Spin Decomposition Theory
!>
!> For two spin-1/2 particles, the Clebsch-Gordan decomposition gives:
!>
!>     ½ ⊗ ½ = 0 ⊕ 1  (singlet ⊕ triplet)
!>
!> The α-β block of the 2-RDM decomposes as: D²_αβ = D²_singlet + D²_triplet
!>
!> For a **true singlet state** (total spin S=0), the triplet component
!> vanishes identically and D²_αβ ≡ D²_singlet.
!>
!> ## Applications
!>
!>   - **N-representability:** Singlet and triplet blocks have independent
!>     positivity constraints (D-condition, Q-condition) in N-representability theory.
!>
!>   - **Purification:** Correction algorithms process singlet/triplet separately
!>     to preserve spin symmetry while enforcing positivity.
!>
!>   - **Symmetry Adaptation:** TD-2RDM propagation can exploit spin symmetry
!>     by evolving only the relevant spin block for singlet states.
!>
!> **Complexity:** O(nOa⁴) with O(nOa⁴) output allocation.
!>
!> @param[out]    rdm2singlet  Singlet-projected 2-RDM (allocated if needed), shape (nOa, nOa, nOa, nOa).
!> @param[in]     rdm2ab       Input α-β 2-RDM block, shape (nOa, nOa, nOa, nOa).
!>                             Convention: rdm2ab(iα, iβ, jα, jβ) with α-orbitals in dims 1,3.
!>
!> @see RdmAlgebra_FillTripletPart for the complementary triplet extraction.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmAlgebra_FillSingletPart(rdm2singlet, rdm2ab)

    complex(R64), intent(out), allocatable :: rdm2singlet(:, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    integer(I32) :: i1, i2, j1, j2, nOa

    nOa = size(rdm2ab, 1)

    if (.not. allocated(rdm2singlet)) allocate (rdm2singlet(nOa, nOa, nOa, nOa))

    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa

            rdm2singlet(i1, i2, j1, j2) = 0.25_R64 * (rdm2ab(i1, i2, j1, j2) + rdm2ab(i2, i1, j2, j1) + &
                                                      rdm2ab(i1, i2, j2, j1) + rdm2ab(i2, i1, j1, j2))

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Extract triplet (S=1) component from α-β 2-RDM.
!>
!> @details
!> Decomposes the spin-resolved α-β 2-RDM into spin irreducible representations
!> (irreps) under the SU(2) spin rotation group. The triplet (S=1) component
!> is **antisymmetric** under combined exchange of spatial and spin indices:
!>
!>     D²_triplet(i1,i2,j1,j2) = (1/4) * [ D²_αβ(i1,i2,j1,j2) + D²_αβ(i2,i1,j2,j1)
!>                                       - D²_αβ(i1,i2,j2,j1) - D²_αβ(i2,i1,j1,j2) ]
!>
!> ## Structure of the Triplet Block
!>
!> The triplet component has special structure arising from Pauli exclusion:
!>
!>   - **Diagonal zeros:** D²_triplet(i,i,j,j) = 0 (same-orbital elements vanish)
!>   - **Off-diagonal antisymmetry:** D²_triplet(i1,i2,j1,j2) = -D²_triplet(i2,i1,j1,j2)
!>
!> This routine exploits these zeros by skipping i1=i2 and j1=j2 loops.
!>
!> ## Applications
!>
!>   - **True Singlet Detection:** For S=0 states, D²_triplet ≡ 0, providing
!>     a numerical check for singlet character: ||D²_triplet|| ≈ 0.
!>
!>   - **Triplet States (S=1, Sz=0):** The α-β block contains nonzero triplet
!>     contribution; this projection isolates that component.
!>
!>   - **Purification:** N-representability corrections must be applied to
!>     singlet and triplet blocks separately to preserve spin symmetry.
!>
!> **Complexity:** O(nOa⁴) with O(nOa⁴) output allocation.
!>               Skips diagonal elements, effective complexity ~ O(nOa²(nOa-1)²).
!>
!> @param[out]    rdm2triplet  Triplet-projected 2-RDM (allocated if needed), shape (nOa, nOa, nOa, nOa).
!>                             Elements with i1=i2 or j1=j2 are zero.
!> @param[in]     rdm2ab       Input α-β 2-RDM block, shape (nOa, nOa, nOa, nOa).
!>                             Convention: rdm2ab(iα, iβ, jα, jβ) with α-orbitals in dims 1,3.
!>
!> @see RdmAlgebra_FillSingletPart for the complementary singlet extraction.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmAlgebra_FillTripletPart(rdm2triplet, rdm2ab)

    complex(R64), intent(out), allocatable :: rdm2triplet(:, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    integer(I32) :: i1, i2, j1, j2, nOa

    nOa = size(rdm2ab, 1)

    if (.not. allocated(rdm2triplet)) allocate (rdm2triplet(nOa, nOa, nOa, nOa))

    rdm2triplet(:, :, :, :) = 0.0_R64

    do j2 = 1, nOa
      do j1 = 1, nOa
        if (j1 .ne. j2) then
          do i2 = 1, nOa
            do i1 = 1, nOa
              if (i1 .ne. i2) then

                rdm2triplet(i1, i2, j1, j2) = 0.25_R64 * (rdm2ab(i1, i2, j1, j2) + rdm2ab(i2, i1, j2, j1) - &
                                                          rdm2ab(i1, i2, j2, j1) - rdm2ab(i2, i1, j1, j2))

              end if
            end do
          end do
        end if
      end do
    end do

  end subroutine

end module
