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

    OrbsInit_Load_dir = Json_Get("orbsInit.load.dir", "")
    OrbsInit_Load_lmaxIn = Json_Get("orbsInit.load.lmaxIn", -1)
    if (0 <= OrbsInit_Load_lmaxIn .and. .not. Json_GetExistence("grid.ylm")) &
      error stop "orbsInit.load.lmaxIn requires grid.ylm"

    OrbsInit_InitializeOrb => InitializeOrb

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Load a single orbital from a binary file.
!>
!> @details
!> Reads orbital data from file "orb{bt}_{ind}.in" (in orbsInit.load.dir if
!> given) using DataStorage routines. The file must contain a complex(R64)
!> array matching the grid size, or, with orbsInit.load.lmaxIn, the orbital on
!> a Ylm grid with the same radial grid and mmax but lmax = lmaxIn: the channel
!> layout is l-major, so the stored orbital is a prefix of the current one and
!> the missing high-l channels are zero.
!>
!> @param[out] orb  Complex orbital vector filled from file.
!> @param[in]  ind  Orbital index (used in filename).
!> @param[in]  bt_  Body type (used in filename; defaults to 1 if absent).
  subroutine InitializeOrb(orb, ind, bt_)
    use M_Utils_DataStorage
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous :: orb(:)
    integer(I32), intent(in)              :: ind
    integer(I32), intent(in), optional    :: bt_

    integer(I32)       :: bt, nIn
    character(len=256) :: filename

    bt = 1
    if (present(bt_)) bt = bt_

    write (filename, '("orb",I2.2,"_",I2.2,".in")') bt, ind
    if (len_trim(OrbsInit_Load_dir) > 0) filename = trim(OrbsInit_Load_dir)//"/"//trim(filename)

    if (OrbsInit_Load_lmaxIn < 0) then
      call LoadData(orb, trim(filename), storage_size(orb))
      return
    end if

    if (OrbsInit_Load_lmaxIn > Grid_Ylm_lmax) error stop "orbsInit.load.lmaxIn exceeds grid.ylm.lmax"
    nIn = Grid_Ylm_NChannels(OrbsInit_Load_lmaxIn, Grid_Ylm_mmax) * Grid_Ylm_nRadial
    orb(:) = 0.0_R64
    call LoadData(orb(1:nIn), trim(filename), storage_size(orb))

  end subroutine

end submodule
