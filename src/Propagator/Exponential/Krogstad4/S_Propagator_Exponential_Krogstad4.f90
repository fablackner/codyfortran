! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Propagator_Exponential_Krogstad4) S_Propagator_Exponential_Krogstad4

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Propagator_Exponential_Krogstad4_Fabricate
    use M_Utils_Say
    use M_Propagator
    use M_Propagator_Exponential

    call Say_Fabricate("propagator.exponential.krogstad4")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Propagator_Propagate => Propagator_Exponential_KrogstadStep

  end subroutine

end submodule
