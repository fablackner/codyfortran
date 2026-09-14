! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Reconstruction_Singlet) S_Reconstruction_Singlet

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Singlet_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Reconstruction_Singlet_Ny
    use M_Reconstruction_Singlet_Valdemoro

    call Say_Fabricate("reconstruction.singlet")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("reconstruction.singlet.singletNy")) then
      call Reconstruction_Singlet_Ny_Fabricate

    else if (Json_GetExistence("reconstruction.singlet.valdemoro")) then
      call Reconstruction_Singlet_Valdemoro_Fabricate

    else
      error stop "reconstruction.singlet is missing one of: singletNy, valdemoro"
    end if

  end subroutine

end submodule

