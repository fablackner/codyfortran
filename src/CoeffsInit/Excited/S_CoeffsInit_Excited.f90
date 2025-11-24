! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_CoeffsInit_Excited) S_CoeffsInit_Excited

  implicit none

  integer(I32), allocatable :: creates(:)
  integer(I32), allocatable :: destroys(:)
  integer(I32) :: bodyType1 = 1
  integer(I32) :: bodyType2 = 2

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine CoeffsInit_Excited_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_CoeffsInit

    call Say_Fabricate("coeffsInit.excited")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    creates = Json_Get("coeffsInit.excited.creates", [1])
    destroys = Json_Get("coeffsInit.excited.destroys", [1])
    bodyType1 = Json_Get("coeffsInit.excited.bodyType1", 1)
    bodyType2 = Json_Get("coeffsInit.excited.bodyType2", 2)

    ! Map the local implementation to the main procedure pointer
    CoeffsInit_Initialize => Initialize

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Initialize(coeffs)
    use M_Coeffs

    complex(R64), intent(out), contiguous :: coeffs(:)

    ! Initialize with ground state
    coeffs(:) = 0.0_R64
    coeffs(1) = 1.0_R64

    ! Apply excitations
    call Coeffs_ApplyExcitation(coeffs, creates, destroys, bodyType1)
    call Coeffs_ApplyExcitation(coeffs, creates, destroys, bodyType2)

  end subroutine

end submodule
