! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Operator application interfaces for grid-based many-body methods.
!>
!> This module declares procedure pointers and abstract interfaces for applying
!> kinetic, potential, and interaction operators directly on grid-based
!> representations of many-body states (as 1D or 2D slices). Implementations
!> are bound at runtime by `Method_Mb_GridBased_Fabricate`.
module M_Method_Mb_GridBased
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds grid-based operator application procedures at runtime based on
    !> the chosen configuration and input parameters.
    module subroutine Method_Mb_GridBased_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the procedure for applying the kinetic energy operator.
  procedure(I_Method_Mb_GridBased_ApplyKineticOp), pointer :: Method_Mb_GridBased_ApplyKineticOp
  abstract interface
    !> Applies the Laplacian differential operator to a 1D slice of the many-body
    !> wavefunction. Typically used to calculate the kinetic energy term.
    subroutine I_Method_Mb_GridBased_ApplyKineticOp(dSlice, slice, time, bt_)
      import :: I32, R64
      !> In/out: derivative slice after applying Laplacian (accumulate).
      complex(R64), intent(inout), contiguous :: dSlice(:)
      !> Input: slice to which the Laplacian is applied.
      complex(R64), intent(in), contiguous :: slice(:)
      !> Current time, allowing for time-dependence.
      real(R64), intent(in) :: time
      !> Optional body-type index for mass scaling or operator selection.
      integer(I32), intent(in), optional :: bt_
    end subroutine
  end interface

  !> Pointer to the procedure for applying the potential energy operator.
  procedure(I_Method_Mb_GridBased_ApplyPotentialOp), pointer :: Method_Mb_GridBased_ApplyPotentialOp
  abstract interface
    !> Applies the potential energy operator to a 1D slice of the many-body
    !> wavefunction and accumulates the result (dSlice += V(slice)).
    subroutine I_Method_Mb_GridBased_ApplyPotentialOp(dSlice, slice, time, bt_)
      import :: I32, R64
      !> In/out: where the result is added. Initialize before calling.
      complex(R64), intent(inout), contiguous :: dSlice(:)
      !> Input: slice to which the potential operator is applied.
      complex(R64), intent(in), contiguous :: slice(:)
      !> Current time, allowing for time-dependence.
      real(R64), intent(in) :: time
      !> Optional body-type index for potential selection.
      integer(I32), intent(in), optional :: bt_
    end subroutine
  end interface

  !> Pointer to the procedure for applying the interaction energy operator.
  procedure(I_Method_Mb_GridBased_ApplyInteractionOp), pointer :: Method_Mb_GridBased_ApplyInteractionOp
  abstract interface
    !> Applies the interaction operator to a 2D slice of the many-body wavefunction
    !> representing the interaction between two particles and accumulates the result.
    subroutine I_Method_Mb_GridBased_ApplyInteractionOp(dSlice2d, slice2d, bt1, bt2, time)
      import :: I32, R64
      !> In/out: 2D slice where the interaction result is added.
      complex(R64), intent(inout), contiguous :: dSlice2d(:, :)
      !> Input: 2D slice for the interaction calculation.
      complex(R64), intent(in), contiguous :: slice2d(:, :)
      !> Body-type index for first particle.
      integer(I32), intent(in) :: bt1
      !> Body-type index for second particle.
      integer(I32), intent(in) :: bt2
      !> Current time, allowing for time-dependence.
      real(R64), intent(in) :: time
    end subroutine
  end interface

end module
