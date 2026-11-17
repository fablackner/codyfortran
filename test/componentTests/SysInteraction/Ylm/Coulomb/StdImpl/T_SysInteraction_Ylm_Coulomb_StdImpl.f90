!> Radial Poisson solve of sysInteraction.ylm.coulomb.stdImpl (Green's
!> function sums over the quadrature points) for the manufactured l = 1 pair
!> of YlmCoulombManufactured.inc; measured average error 7.8e-3.
program T_SysInteraction_Ylm_Coulomb_StdImpl
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/componentTests/SysInteraction/Ylm/Coulomb/StdImpl/T_SysInteraction_Ylm_Coulomb_StdImpl.json"

  call Say_Hello
  call CheckManufacturedPoisson(jsonFileName, 1e-2_R64)
  call Say_Goodbye

contains

  include '../YlmCoulombManufactured.inc'

end program
