! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Triples amplitudes: source / D with the orbital-energy denominator D (second order).
!>
!> **JSON path**: `reconstruction.singlet.triples.moellerPlesset`
module M_Reconstruction_Singlet_Triples_MoellerPlesset
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds Reconstruction_Singlet_Triples_FillAmplitudes.
    module subroutine Reconstruction_Singlet_Triples_MoellerPlesset_Fabricate
    end subroutine
  end interface

end module
