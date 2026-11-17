!> Radial Poisson solve of sysInteraction.ylm.coulomb.fullEq (FEDVR equation,
!> LU per call) for the manufactured l = 1 pair of YlmCoulombManufactured.inc;
!> measured average error 4.7e-8.
program T_SysInteraction_Ylm_Coulomb_FullEq
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/componentTests/SysInteraction/Ylm/Coulomb/FullEq/T_SysInteraction_Ylm_Coulomb_FullEq.json"

  call Say_Hello
  call CheckManufacturedPoisson(jsonFileName, 1e-7_R64)
  call Say_Goodbye

contains

  include '../YlmCoulombManufactured.inc'

end program
