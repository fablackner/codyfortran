! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Implementation of the eigen-expansion propagator backend.
!>
!> Propagates quantum states by expanding in the eigenbasis of Ĥ:
!>   |Ψ(t₁)⟩ = Σᵢ cᵢ exp(−iEᵢΔt) |φᵢ⟩
!>
!> where cᵢ = ⟨φᵢ|Ψ(t₀)⟩ are the expansion coefficients.
!>
!> Best suited for:
!> - Time-independent Hamiltonians where diagonalization cost is amortized
!> - Small to medium Hilbert spaces (dim ≲ 10⁴)
!> - High-accuracy requirements over long time intervals
!>
!> Limitations:
!> - Incomplete eigenbasis: If nFound < dim, the projection is approximate
!> - Memory: Stores full eigenvector matrix O(dim × nFound)
!> - Setup cost: Single O(dim³) diagonalization at initialization
!>
!> Requirements:
!> - DiagonalizerList(1) must be fabricated with ApplyMatOnVec bound to Ĥ
submodule(M_Propagator_EigenExpansion) S_Propagator_EigenExpansion

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Binds the Propagator facade pointers to this backend.
!>
!> After calling this routine:
!> - `Propagator_Setup` → Setup (diagonalizes Ĥ once)
!> - `Propagator_Propagate` → Propagate (applies phase factors)
  module subroutine Propagator_EigenExpansion_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Propagator

    call Say_Fabricate("propagator.eigenExpansion")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Propagator_Propagate => Propagate
    Propagator_Setup => Setup

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Diagonalizes the Hamiltonian and stores eigenpairs for later propagation.
!>
!> Calls DiagonalizerList(1)%Diagonalize at t=0 with eigenvector storage.
!> The eigenvalues Eᵢ and eigenvectors |φᵢ⟩ are cached in DiagonalizerList.
  subroutine Setup()
    use M_Utils_Say
    use M_DiagonalizerList

    call Say_Setup("propagator.eigenExpansion")

    ! Diagonalize Ĥ at t=0; evecsQ=.true. requests eigenvector storage
    call DiagonalizerList(1) % e % Diagonalize(0.0_R64, .true.)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Propagates the state using the cached eigenbasis.
!>
!> Algorithm:
!>   1. Project: cᵢ = ⟨φᵢ|Ψ(t₀)⟩ for i = 1..nFound
!>   2. Evolve:  |Ψ(t₁)⟩ = Σᵢ cᵢ exp(−iEᵢΔt) |φᵢ⟩
!>
!> Complexity: O(nFound × dim) per propagation step
!>
!> Note: Uses only the eigenpairs found by the diagonalizer. If nFound < dim,
!> components orthogonal to the stored eigenbasis are projected out.
!>
!> Arguments:
!>   state  - Complex state vector, modified in-place
!>   t0     - Initial time (used only for Δt = t1 − t0)
!>   t1     - Final time
  subroutine Propagate(state, t0, t1)
    use M_Utils_Constants
    use M_DiagonalizerList

    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in)    :: t0
    real(R64), intent(in)    :: t1

    integer(I32) :: i, nEvals
    real(R64) :: dt
    complex(R64), allocatable :: c(:)

    dt = t1 - t0
    nEvals = DiagonalizerList(1) % e % nFound

    allocate (c(nEvals))

    ! Project state onto eigenbasis: cᵢ = ⟨φᵢ|Ψ⟩
    do i = 1, nEvals
      c(i) = dot_product(DiagonalizerList(1) % e % evecs(:, i), state)
    end do

    ! Reconstruct with phase factors: |Ψ(t₁)⟩ = Σᵢ cᵢ exp(−iEᵢΔt) |φᵢ⟩
    state(:) = 0.0_R64
    do i = 1, nEvals
      associate (diagonalizer => DiagonalizerList(1) % e)
        state(:) = state(:) + exp(-IU * diagonalizer % evals(i) * dt) * c(i) * diagonalizer % evecs(:, i)
      end associate
    end do

  end subroutine

end submodule
