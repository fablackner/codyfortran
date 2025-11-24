! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Orbs) S_Orbs

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Orbs_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Method_Mb_OrbBased

    call Say_Fabricate("orbs")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Orbs_ProjectOnSubspace => ProjectOnSubspace
    Orbs_Orthonormalize => Orthonormalize
    Orbs_SaveOrbs => SaveOrbs

    Orbs_restrictedQ = Json_Get("orbs.restrictedQ", .false.)

    if (Orbs_restrictedQ) Orbs_nOrbsInState = Method_Mb_OrbBased_nOrbsSum / 2
    if (.not. Orbs_restrictedQ) Orbs_nOrbsInState = Method_Mb_OrbBased_nOrbsSum

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Orthonormalize(orbs)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Grid

    complex(R64), intent(inout), contiguous  :: orbs(:, :)

    integer(I32) :: ibt, startOrb, endOrb

    startOrb = 1
    do ibt = 1, Method_Mb_nBodyTypes

      endOrb = startOrb + Method_Mb_OrbBased_nOrbs(ibt) - 1

      call Grid_Orthonormalize(orbs(:, startOrb:endOrb))

      startOrb = endOrb + 1

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ProjectOnSubspace(dOrbs, orbs)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Grid

    complex(R64), intent(inout), contiguous  :: dOrbs(:, :)
    complex(R64), intent(in), contiguous     :: orbs(:, :)

    integer(I32) :: ibt, startOrb, endOrb

    startOrb = 1
    do ibt = 1, Method_Mb_nBodyTypes

      endOrb = startOrb + Method_Mb_OrbBased_nOrbs(ibt) - 1

      call Grid_ProjectOnSubspace(dOrbs(:, startOrb:endOrb), orbs(:, startOrb:endOrb))

      startOrb = endOrb + 1

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine SaveOrbs(orbs)
    use M_Utils_DataStorage
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(in), contiguous  :: orbs(:, :)

    integer(I32) :: index, ibt, i1
    character(len=256) :: filename

    i1 = 0
    do ibt = 1, Method_Mb_nBodyTypes
      do index = 1, Method_Mb_OrbBased_nOrbs(ibt)
        i1 = i1 + 1

        write (filename, '("orb",I2.2,"_",I2.2,".dat")') ibt, index
        call SaveData(orbs(:, i1), trim(filename), storage_size(orbs(:, i1)))

      end do
    end do

  end subroutine

end submodule
