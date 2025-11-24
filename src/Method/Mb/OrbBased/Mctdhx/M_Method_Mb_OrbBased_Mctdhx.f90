! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> MCTDHX, orbital-based method with time-dependent orbitals and coefficients.
!>
!> Both CI coefficients and single-particle orbitals evolve in time according to
!> the Dirac–Frenkel time-dependent variational principle, enabling a compact yet
!> accurate representation of the many-body wavefunction compared to static-orbital
!> approaches. Procedure pointers are bound by the fabricate routine.
module M_Method_Mb_OrbBased_Mctdhx
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds MCTDHX procedure pointers at runtime and initializes state layout
    !> for dynamic orbitals and coefficients. Typically sets:
    !> - Method_Setup, Method_TimeDerivative
    !> - Method_Mb_OrbBased_TimeDerivativeOrbsLin
    !> - Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
    module subroutine Method_Mb_OrbBased_Mctdhx_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the procedure for extracting CI coefficients from the packed state vector.
  procedure(I_Method_Mb_OrbBased_GetCoeffs), pointer :: M_Method_Mb_OrbBased_Mctdhx_GetCoeffs
  abstract interface
    !> Returns a pointer view onto the CI coefficients contained in the method
    !> state vector. No allocation or copy is performed.
    !>
    !> Lifetime
    !> - The returned pointer becomes invalid if the underlying `state` is
    !>   reallocated or its layout changes during subsequent setup calls.
    subroutine I_Method_Mb_OrbBased_GetCoeffs(coeffs, state)
      import :: R64
      complex(R64), pointer, intent(out) :: coeffs(:)
      complex(R64), intent(in), contiguous, target :: state(:)
    end subroutine
  end interface

end module
