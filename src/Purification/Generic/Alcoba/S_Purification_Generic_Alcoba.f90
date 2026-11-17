! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Purification_Generic_Alcoba) S_Purification_Generic_Alcoba

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Purification_Generic_Alcoba_Fabricate
    use M_Utils_Say

    call Say_Fabricate("purification.generic.alcoba")

    ! Placeholder: without this stop the run would continue unpurified
    error stop "purification.generic.alcoba is not implemented"

  end subroutine

end submodule

