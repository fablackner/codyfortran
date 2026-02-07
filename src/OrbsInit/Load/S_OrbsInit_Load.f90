! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_OrbsInit_Load) S_OrbsInit_Load

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine OrbsInit_Load_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit

    call Say_Fabricate("orbsInit.load")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_InitializeOrb => InitializeOrb

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine InitializeOrb(orb, ind, bt_)
    use M_Utils_DataStorage

    complex(R64), intent(out), contiguous        :: orb(:)
    integer(I32), intent(in)         :: ind
    integer(I32), intent(in), optional         :: bt_

    integer(I32)        :: bt
    character(len=256) :: filename

    if (.not. present(bt_)) bt = 1
    if (present(bt_)) bt = bt_

    write (filename, '("orb",I2.2,"_",I2.2,".in")') bt_, ind

    call LoadData(orb, trim(filename), storage_size(orb))

  end subroutine

end submodule
