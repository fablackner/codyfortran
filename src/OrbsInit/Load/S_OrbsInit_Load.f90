! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for file-based orbital loading.
submodule(M_OrbsInit_Load) S_OrbsInit_Load

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Bind the InitializeOrb pointer for file-based loading.
  module subroutine OrbsInit_Load_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit

    call Say_Fabricate("orbsInit.load")

    OrbsInit_InitializeOrb => InitializeOrb

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Load a single orbital from a binary file.
!>
!> @details
!> Reads orbital data from file "orb{bt}_{ind}.in" using DataStorage routines.
!> The file must contain a complex(R64) array matching the grid size.
!>
!> @param[out] orb  Complex orbital vector filled from file.
!> @param[in]  ind  Orbital index (used in filename).
!> @param[in]  bt_  Body type (used in filename; defaults to 1 if absent).
  subroutine InitializeOrb(orb, ind, bt_)
    use M_Utils_DataStorage

    complex(R64), intent(out), contiguous :: orb(:)
    integer(I32), intent(in)              :: ind
    integer(I32), intent(in), optional    :: bt_

    integer(I32)       :: bt
    character(len=256) :: filename

    bt = 1
    if (present(bt_)) bt = bt_

    write (filename, '("orb",I2.2,"_",I2.2,".in")') bt, ind

    call LoadData(orb, trim(filename), storage_size(orb))

  end subroutine

end submodule
