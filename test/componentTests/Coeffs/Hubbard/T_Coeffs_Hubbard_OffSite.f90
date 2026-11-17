!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> coeffs.hubbard refuses a two-body Hamiltonian it cannot represent: h2 of
!> the interaction V n_{1,up} n_{2,dn} between neighbouring sites.
!> ctest-pass-regex: h2 is not an on-site interaction
program T_Coeffs_Hubbard_OffSite
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Method
  use M_Method_Mb_OrbBased
  use M_Coeffs

  implicit none

  complex(R64), allocatable :: h2(:, :, :, :), c(:), d(:)
  integer(I32) :: nO, oUp, oDn

  call Say_Hello

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Coeffs/Hubbard/T_Coeffs_Hubbard_OffSite.json", &
                         relativeToProjectDirQ_=.true.)
  call Grid_Fabricate
  call Method_Fabricate
  call Coeffs_Fabricate
  call Grid_Setup
  call Coeffs_Setup

  nO = Method_Mb_OrbBased_nOrbsSum
  oUp = Method_Mb_OrbBased_nOrbsStart(1) - 1
  oDn = Method_Mb_OrbBased_nOrbsStart(2) - 1
  allocate (h2(nO, nO, nO, nO))
  h2(:, :, :, :) = 0.0_R64
  h2(oUp + 1, oDn + 2, oUp + 1, oDn + 2) = 0.5_R64
  h2(oDn + 2, oUp + 1, oDn + 2, oUp + 1) = 0.5_R64

  allocate (c(Coeffs_nCoeffs), d(Coeffs_nCoeffs))
  c(:) = 1.0_R64
  d(:) = 0.0_R64
  call Coeffs_ApplyH2FillRdm2(c, dCoeffs_=d, h2_=h2)

  error stop "T_Coeffs_Hubbard_OffSite failure: the off-site interaction was accepted"

end program
