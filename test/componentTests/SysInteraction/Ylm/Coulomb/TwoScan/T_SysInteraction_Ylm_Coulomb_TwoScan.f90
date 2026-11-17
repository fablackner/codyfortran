!> Radial Poisson solve of sysInteraction.ylm.coulomb.twoScan (the sums of
!> stdImpl in two prefix scans) for the manufactured l = 1 pair of
!> YlmCoulombManufactured.inc; measured average error 7.8e-3, as stdImpl.
program T_SysInteraction_Ylm_Coulomb_TwoScan
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/componentTests/SysInteraction/Ylm/Coulomb/TwoScan/T_SysInteraction_Ylm_Coulomb_TwoScan.json"

  call Say_Hello
  call CheckManufacturedPoisson(jsonFileName, 1e-2_R64)
  call Say_Goodbye

contains

  include '../YlmCoulombManufactured.inc'

end program
