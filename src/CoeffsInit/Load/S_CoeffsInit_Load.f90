! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_CoeffsInit_Load) S_CoeffsInit_Load

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine CoeffsInit_Load_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_CoeffsInit

    call Say_Fabricate("coeffsInit.load")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Map the local implementation to the main procedure pointer
    CoeffsInit_Initialize => Initialize

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Initialize(coeffs)
    use M_Utils_DataStorage

    complex(R64), intent(out), contiguous :: coeffs(:)
    character(len=256) :: filename

    write (filename, '("coeffs.dat")')

    call LoadData(coeffs, filename, storage_size(coeffs))

  end subroutine

end submodule
