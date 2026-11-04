! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Sb) S_Method_Sb
  !-----------------------------------------------------------------------------
  ! Single-body (Sb) method implementation.
  !
  ! Evolves a single wavefunction ψ(r,t) under the kinetic and external
  ! potential operators and, if a sysGauge block exists, the laser coupling
  ! (SysGauge: velocity gauge on Ylm, length gauge sysGauge.fem.length on
  ! grid.fem). No mean-field, exchange, or correlation terms are included.
  !
  ! The TDSE solved is:
  !   i ∂ψ/∂t = (T̂ + V̂_ext + Ĥ_A(t)) ψ
  !
  ! Useful for:
  !   - Testing grid implementations and discretization schemes
  !   - Studying single-particle phenomena (tunneling, wave packets, scattering)
  !   - Validating propagator accuracy against analytical solutions
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Gauge (laser) term of H ψ, bound in Method_Sb_Fabricate as
  !> Method_Mb_OrbBased_ApplyGaugeOp for the orbital methods: AddGaugeTerm if
  !> a sysGauge block exists, SkipGaugeTerm otherwise.
  procedure(I_AccumulateGaugeTerm), pointer :: AccumulateGaugeTerm => null()
  abstract interface
    !> hPsi += Ĥ_A(time) ψ.
    subroutine I_AccumulateGaugeTerm(hPsi, psi, time)
      import :: R64
      complex(R64), intent(inout), contiguous :: hPsi(:)
      complex(R64), intent(in), contiguous    :: psi(:)
      real(R64), intent(in)                   :: time
    end subroutine
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Sb_Fabricate
    !---------------------------------------------------------------------------
    ! Binds single-body procedure pointers.
    !
    ! After this call:
    !   - Method_Setup          → Setup (allocates state, initializes orbital)
    !   - Method_GetEnergy      → GetEnergy (computes ⟨ψ|H|ψ⟩)
    !   - Method_TimeDerivative → TimeDerivative (computes -i·H·ψ)
    !   - the gauge term of H   → AddGaugeTerm if "sysGauge" exists
    !---------------------------------------------------------------------------
    use M_Utils_Json
    use M_Utils_Say
    use M_Method

    call Say_Fabricate("method.sb")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Method_Setup => Setup
    Method_GetEnergy => GetEnergy
    Method_TimeDerivative => TimeDerivative

    ! Laser coupling only if configured (SysGauge_Fabricate binds the operator)
    if (Json_GetExistence("sysGauge")) then
      AccumulateGaugeTerm => AddGaugeTerm
    else
      AccumulateGaugeTerm => SkipGaugeTerm
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    !---------------------------------------------------------------------------
    ! Allocates the single-body state vector and initializes the wavefunction.
    !
    ! State layout: Method_state(1:nPoints) = ψ(r) on the spatial grid.
    ! Initialization is delegated to OrbsInit_InitializeOrb with orbital index 1.
    !---------------------------------------------------------------------------
    use M_Utils_Say
    use M_Grid
    use M_OrbsInit
    use M_Method

    call Say_Setup("method.sb")

    allocate (Method_state(1:Grid_nPoints))

    ! Initialize the single-body wavefunction from OrbsInit configuration
    call OrbsInit_InitializeOrb(Method_state, 1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function GetEnergy(time) result(res)
    !---------------------------------------------------------------------------
    ! Computes the energy expectation value E = ⟨ψ|Ĥ|ψ⟩.
    !
    ! Uses the identity H|ψ⟩ = i·(d|ψ⟩/dt), so:
    !   E = Re[⟨ψ| i·(dψ/dt)⟩] = Re[⟨ψ| H|ψ⟩]
    !---------------------------------------------------------------------------
    use M_Utils_Constants
    use M_Grid
    use M_Method

    real(R64)                :: res
    real(R64), intent(in)    :: time

    complex(R64), allocatable :: Hpsi(:)

    ! Compute H|ψ⟩ via the time derivative: dψ/dt = -i·H·ψ  →  H·ψ = i·dψ/dt
    allocate (Hpsi(Grid_nPoints))

    call TimeDerivative(Hpsi, Method_state, time)
    Hpsi = IU * Hpsi  ! H|ψ⟩ = i·(dψ/dt)

    ! Energy = Re[⟨ψ|H|ψ⟩] (imaginary part vanishes for Hermitian H)
    res = real(Grid_InnerProduct(Method_state, Hpsi), kind=R64)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivative(dState, state, time)
    !---------------------------------------------------------------------------
    ! Computes the right-hand side of the single-body TDSE:
    !   dψ/dt = -i·(T̂ + V̂_ext + Ĥ_A(t))·ψ
    !
    ! The Hamiltonian action is computed as:
    !   1. Apply kinetic operator:    dState = T̂·ψ
    !   2. Apply external potential:  dState += V̂_ext·ψ
    !   3. Apply gauge coupling:      dState += Ĥ_A(t)·ψ  (if sysGauge exists)
    !   4. Multiply by -i for TDSE:   dState *= -i
    !---------------------------------------------------------------------------
    use M_Utils_Constants
    use M_SysKinetic
    use M_SysPotential

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), allocatable :: externalPotential(:)
    complex(R64), allocatable :: Vpsi(:)

    ! Initialize accumulator for H·ψ
    dState(:) = 0.0_R64

    allocate (Vpsi, mold=state)

    ! Apply kinetic operator: dState += T̂·ψ
    call SysKinetic_MultiplyWithKineticOp(dState, state, time)

    ! Apply external potential: dState += V̂·ψ
    call SysPotential_FillExternalPotential(externalPotential, time)
    call SysPotential_MultiplyWithExternalPotential(Vpsi, externalPotential, state)
    dState(:) = dState(:) + Vpsi(:)

    ! Apply the laser coupling: dState += Ĥ_A(t)·ψ (no-op without sysGauge)
    call AccumulateGaugeTerm(dState, state, time)

    ! Convert to time derivative: dψ/dt = -i·H·ψ
    dState(:) = -IU * dState(:)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AddGaugeTerm(hPsi, psi, time)
    !---------------------------------------------------------------------------
    ! hPsi += Ĥ_A(time) ψ with the SysGauge operator of body type 1.
    !---------------------------------------------------------------------------
    use M_SysGauge

    complex(R64), intent(inout), contiguous :: hPsi(:)
    complex(R64), intent(in), contiguous    :: psi(:)
    real(R64), intent(in)                   :: time

    complex(R64), allocatable :: gaugePsi(:)

    allocate (gaugePsi, mold=psi)
    call SysGauge_MultiplyWithGaugeOp(gaugePsi, psi, time, bt_=1)
    hPsi(:) = hPsi(:) + gaugePsi(:)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine SkipGaugeTerm(hPsi, psi, time)
    !---------------------------------------------------------------------------
    ! No laser coupling configured.
    !---------------------------------------------------------------------------
    use M_Utils_UnusedVariables

    complex(R64), intent(inout), contiguous :: hPsi(:)
    complex(R64), intent(in), contiguous    :: psi(:)
    real(R64), intent(in)                   :: time

    if (.false.) call UnusedVariables_Mark(hPsi, psi, time)

  end subroutine

end submodule
