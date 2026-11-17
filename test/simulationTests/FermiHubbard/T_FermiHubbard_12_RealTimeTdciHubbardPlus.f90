!> Real-time TDCI of a 5-site Fermi-Hubbard chain with coeffs.hubbard.plusSpinSym
!> (FermiHubbardRealTime.inc).
program T_FermiHubbard_12_RealTimeTdciHubbardPlus
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = "test/simulationTests/FermiHubbard/T_FermiHubbard_12_RealTimeTdciHubbardPlus.json"

  call Say_Hello
  call CheckRealTimeTdci(jsonFileName)
  call Say_Goodbye

contains

  include 'FermiHubbardRealTime.inc'

end program
