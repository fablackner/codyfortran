! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Singlet (restricted) fermionic backend for 2‑RDM operations.
!>
!> This backend exploits spin‑singlet symmetry in a spatial‑orbital description.
!> It provides conversion helpers between the backend’s linearized 2‑RDM and an
!> alpha–beta (ab) 4‑index tensor over spatial orbitals. Concrete wiring is done by
!> `TwoRdm_Fermionic_Singlet_Fabricate`.
module M_TwoRdm_Fermionic_Singlet
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Activate the singlet fermionic backend and wire procedure pointers.
    !>
    !> Responsibilities:
    !> - Assign `M_TwoRdm` procedure pointers to singlet‑optimized implementations.
    !> - Define the linearization layout consistent with the ab spatial tensor mapping.
    module subroutine TwoRdm_Fermionic_Singlet_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! local data
  !=============================================================================

  !> Pointer to the procedure converting a 4D ab spatial‑orbital 2‑RDM to linearized storage.
  procedure(I_TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab), pointer :: TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab
  abstract interface
    !> Pack the ab spatial‑orbital 2‑RDM tensor into the backend’s linearized format.
    !>
    !> The ab tensor is a 4‑index object over spatial orbitals; spin structure is
    !> implicit in the singlet mapping adopted by this backend.
    subroutine I_TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab(twoRdmState, rdm2ab)
      import :: R64
      !> Output linearized 2‑RDM (pre‑allocated; size must match backend layout).
      complex(R64), intent(out), contiguous  :: twoRdmState(:)
      !> Input ab 2‑RDM tensor with shape (n, n, n, n) over spatial orbitals.
      complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)
    end subroutine
  end interface

  !> Pointer to the procedure converting the linearized 2‑RDM to a 4D ab spatial tensor.
  procedure(I_TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState), pointer :: TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState
  abstract interface
    !> Unpack the backend’s linearized 2‑RDM into an ab 4‑index spatial‑orbital tensor.
    subroutine I_TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
      import :: R64
      !> Allocated on output: ab tensor with shape (n, n, n, n) over spatial orbitals.
      complex(R64), intent(out), allocatable :: rdm2ab(:, :, :, :)
      !> Input linearized 2‑RDM.
      complex(R64), intent(in), contiguous   :: twoRdmState(:)
    end subroutine
  end interface

end module
