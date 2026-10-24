!> Method_Mb_OrbBased_OrthonormalizeKeepingState recovers the orthonormal
!> orbitals and CI vector of a state rewritten with non-orthonormal orbitals
!> (He1dOrthonormalize.inc); restricted orbitals.
program T_He1d_10_OrthonormalizeKeepingStateRestricted
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  call Say_Hello
  call CheckOrthonormalizeKeepingState("test/simulationTests/He1d/T_He1d_10_OrthonormalizeKeepingStateRestricted.json")
  call Say_Goodbye

contains

  include 'He1dOrthonormalize.inc'

end program
