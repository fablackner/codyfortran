! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Utilities to apply the problem Hamiltonian and time-derivative operators.
!>
!> This module provides thin wrappers that delegate to the active Method
!> implementation (via `M_Method`) and arrange the sign/convention so that
!> callers can either obtain the time derivative d|ψ>/dt or the Hamiltonian
!> action H|ψ> as needed.
!>
!> Notes
!> - `Actions_TimeDerivative` computes d|ψ>/dt as provided by the selected method.
!> - `Actions_ImagTimeDerivative` computes the imaginary-time derivative.
!> - `Actions_ApplyHamiltonian` returns H|ψ> by converting the time derivative
!>   with the Schrödinger convention d|ψ>/dt = -i H|ψ> (i.e., H|ψ> = i d|ψ>/dt).
module M_Utils_Actions
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Apply the Hamiltonian to a state vector.
  !>
  !> This routine calls the active method to obtain the time derivative
  !> dState := d|ψ>/dt, then converts it to the Hamiltonian action using
  !> H|ψ> = i d|ψ>/dt. The output `dState` is thus the Hamiltonian applied
  !> to the input `state` at the given `time`.
  !>
  !> Parameters
  !> - dState [out]: complex vector receiving H|ψ> (same shape as `state`).
  !> - state  [in]:  complex state vector |ψ>.
  !> - time   [in]:  physical time at which to evaluate the derivative.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Actions_ApplyHamiltonian(dState, state, time)
    use M_Utils_Constants
    use M_Method

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in) :: time

    call Method_TimeDerivative(dState, state, time)

    dState(:) = IU * dState(:) ! out is the Hamiltonian applied to the state

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Actions_TimeDerivative(dState, state, time)
    use M_Utils_Constants
    use M_Method

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in) :: time

    ! This is a wrapper for the time derivative
    call Method_TimeDerivative(dState, state, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Actions_ImagTimeDerivative(dState, state, time)
    use M_Utils_Constants
    use M_Method

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in) :: time

    call Method_TimeDerivative(dState, state, time)

    dState(:) = -IU * dState(:) ! out is the imaginary time derivative

  end subroutine

end module
