! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Reconstruction_SingletTr0) S_Reconstruction_SingletTr0

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_SingletTr0_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Reconstruction_SingletTr0_Valdemoro
    use M_Reconstruction_SingletTr0_Ny

    call Say_Fabricate("reconstruction.singletTr0")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("reconstruction.singletTr0.valdemoro")) then
      call Reconstruction_SingletTr0_Valdemoro_Fabricate

    else if (Json_GetExistence("reconstruction.singletTr0.ny")) then
      call Reconstruction_SingletTr0_Ny_Fabricate

    else
      error stop "reconstruction.singletTr0 is missing one of: valdemoro, ny"
    end if

  end subroutine

end submodule
