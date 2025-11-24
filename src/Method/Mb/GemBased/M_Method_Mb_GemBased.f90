! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Operator application interfaces for GEM-based many-body methods.
!>
!> This module declares procedure pointers and abstract interfaces for applying
!> one- and two-body operators in a basis represented by GEM (Gaussian or other
!> basis-expansion) functions. Concrete implementations are bound at runtime by
!> `Method_Mb_GemBased_Fabricate` depending on the selected configuration.
module M_Method_Mb_GemBased
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds GEM-based operator application procedures at runtime based on input
    !> parameters and selected configuration.
    module subroutine Method_Mb_GemBased_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the procedure for applying the kinetic energy operator.
  procedure(I_Method_Mb_GemBased_ApplyKineticOp), pointer :: Method_Mb_GemBased_ApplyKineticOp
  abstract interface
    !> Applies the Laplacian (kinetic) operator to a set of orbitals in the GEM basis.
    !>
    !> Conventions
    !> - Arrays are column-major; columns of `orbs` represent single-particle orbitals.
    !> - The contribution is accumulated in `dOrbs` (i.e., dOrbs += K(orbs)).
    subroutine I_Method_Mb_GemBased_ApplyKineticOp(dOrbs, orbs, time)
      import :: R64
      !> In/out: derivative orbitals after applying Laplacian (accumulate).
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      !> Input: orbitals to which the Laplacian is applied.
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Current time, allowing for time-dependence.
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the procedure for applying the potential energy operator.
  procedure(I_Method_Mb_GemBased_ApplyPotentialOp), pointer :: Method_Mb_GemBased_ApplyPotentialOp
  abstract interface
    !> Applies the one-body potential operator to the orbitals and accumulates
    !> the result in `dOrbs` (i.e., dOrbs += V(orbs)).
    subroutine I_Method_Mb_GemBased_ApplyPotentialOp(dOrbs, orbs, time)
      import :: R64
      !> In/out: where the result is accumulated. Initialize before calling.
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      !> Input: orbitals to which the potential operator is applied.
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Current time, allowing for time-dependence.
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the procedure for applying the two-body interaction operator.
  procedure(I_Method_Mb_GemBased_ApplyInteractionOp), pointer :: Method_Mb_GemBased_ApplyInteractionOp
  abstract interface
    !> Applies the two-body interaction operator to the orbitals (e.g., Coulomb
    !> or contact interactions) and accumulates the result in `dOrbs`.
    !>
    !> Notes
    !> - The integers `i2` and `j2` select a method-specific interaction channel
    !>   or tensor block (e.g., body-type pair or spin channel).
    subroutine I_Method_Mb_GemBased_ApplyInteractionOp(dOrbs, orbs, i2, j2, time)
      import :: I32, R64
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      complex(R64), intent(in), contiguous   :: orbs(:, :)
      integer(I32), intent(in)               :: i2
      integer(I32), intent(in)               :: j2
      !> Current time, allowing for time-dependence.
      real(R64), intent(in) :: time
    end subroutine
  end interface

end module
