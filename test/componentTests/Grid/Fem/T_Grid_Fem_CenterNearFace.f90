program T_Grid_Fem_CenterNearFace
  !-----------------------------------------------------------------------------
  ! A center closer than 0.5 hMin to a box face (center 2 lies 0.05 below the
  ! upper z face, hMin 0.2) is rejected at fabrication: the graded mesh would
  ! get a sliver element (or the nucleus) at the Dirichlet surface.
  !-----------------------------------------------------------------------------
!> ctest-pass-regex: center 2 is closer than 0.5.hMin to a box face along z
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid

  implicit none

  call Say_Hello

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Grid/Fem/T_Grid_Fem_CenterNearFace.json", &
                         relativeToProjectDirQ_=.true.)

  ! must stop with the expected message
  call Grid_Fabricate

  print *, "T_Grid_Fem_CenterNearFace: the center next to the box face was not rejected"

end program
