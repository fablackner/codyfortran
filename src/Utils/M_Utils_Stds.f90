! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default no-op callbacks used as safe placeholders.
!>
!> Provides stub implementations for optional hooks so that modules can
!> install behavior selectively without adding conditional logic at call sites.
module M_Utils_NoOpProcedures
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine NoOpProcedures_Setup

    ! this is just a dummy routine that does nothing

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine NoOpProcedures_ApplyAbsorber(orbs)
    use M_Utils_UnusedVariables

    complex(R64), intent(inout), contiguous  :: orbs(:, :)

    ! this is just a dummy routine that does nothing

    if (.false.) call UnusedVariables_Mark(orbs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine NoOpProcedures_ApplyInteractionOp(dOrbs, orbs, i2, j2, time)
    use M_Utils_UnusedVariables
    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    integer(I32), intent(in)                :: i2
    integer(I32), intent(in)                :: j2
    real(R64), intent(in)                   :: time

    ! do nothing, just a dummy routine

    if (.false.) call UnusedVariables_Mark(dOrbs, orbs, i2, j2, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine NoOpProcedures_ApplyPotentialOp(dOrbs, orbs, time)
    use M_Utils_UnusedVariables
    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    real(R64), intent(in) :: time

    ! do nothing, just a dummy routine

    if (.false.) call UnusedVariables_Mark(dOrbs, orbs, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine NoOpProcedures_PurifyTwoRdmState(twoRdmState)
    use M_Utils_UnusedVariables

    complex(R64), intent(inout), contiguous  :: twoRdmState(:)

    ! do NoPurification!

    if (.false.) call UnusedVariables_Mark(twoRdmState)

  end subroutine

end module
