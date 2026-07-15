program T_Grid_Prolate
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Prolate
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: iGrid
  real(R64) :: volume, xi, eta, r1
  complex(R64) :: norm, ovlp
  complex(R64), allocatable :: orb1s(:), orbSet(:, :), src(:)
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/Grid/Prolate/T_Grid_Prolate.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Grid_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  !------------------------------------
  ! metric: integral of (ximax - xi), which vanishes at the excluded Dirichlet
  ! endpoint ximax, against dV = a^3 (xi^2 - eta^2) dxi deta dphi:
  !   I = 2 pi a^3 [ 2X(X^3-1)/3 - (2/3)X(X-1) - (X^4-1)/2 + (X^2-1)/3 ],  X = ximax
  ! The quadrature is exact for this polynomial, so this pins weights and
  ! Jacobian including the endpoint-exclusion bookkeeping.
  !------------------------------------

  associate (x => Grid_Prolate_ximax)
    volume = TWOPI * Grid_Prolate_a**3 * (2.0_R64 * x * (x**3 - 1.0_R64) / 3.0_R64 &
                                          - 2.0_R64 / 3.0_R64 * x * (x - 1.0_R64) &
                                          - (x**4 - 1.0_R64) / 2.0_R64 &
                                          + (x**2 - 1.0_R64) / 3.0_R64)
  end associate

  call check(error, &
             sum(Grid_Prolate_spatialWeights * (Grid_Prolate_ximax - Grid_Prolate_xiCoord(1:Grid_Prolate_nSpatial))) &
             * TWOPI, &
             volume, thr=1e-10_R64 * volume)
  if (allocated(error)) error stop "T_Grid_Prolate metric integral failure"

  !------------------------------------
  ! norm of a hydrogen 1s orbital centered on the focus at z = -a:
  ! psi = exp(-r1)/sqrt(pi) with r1 = a (xi + eta), stored as the m = 0
  ! channel f_0 = sqrt(2 pi) psi
  !------------------------------------

  allocate (orb1s(Grid_nPoints))
  orb1s = (0.0_R64, 0.0_R64)
  do iGrid = 1, Grid_nPoints
    if (Grid_Prolate_mCoord(iGrid) .ne. 0) cycle
    xi = Grid_Prolate_xiCoord(iGrid)
    eta = Grid_Prolate_etaCoord(iGrid)
    r1 = Grid_Prolate_a * (xi + eta)
    orb1s(iGrid) = sqrt(TWOPI) * exp(-r1) / sqrt(PI)
  end do

  norm = Grid_InnerProduct(orb1s, orb1s)
  call check(error, real(norm, kind=R64), 1.0_R64, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Grid_Prolate 1s norm failure"

  !------------------------------------
  ! orthonormalization via the shared Gram-Schmidt
  !------------------------------------

  allocate (orbSet(Grid_nPoints, 2))
  orbSet(:, 1) = orb1s
  orbSet(:, 2) = orb1s * Grid_Prolate_xiCoord(:)  ! linearly independent partner

  call Grid_Orthonormalize(orbSet)

  ovlp = Grid_InnerProduct(orbSet(:, 1), orbSet(:, 2))
  call check(error, abs(ovlp), 0.0_R64, thr=1e-12_R64)
  if (allocated(error)) error stop "T_Grid_Prolate orthonormalize overlap failure"

  norm = Grid_InnerProduct(orbSet(:, 2), orbSet(:, 2))
  call check(error, real(norm, kind=R64), 1.0_R64, thr=1e-12_R64)
  if (allocated(error)) error stop "T_Grid_Prolate orthonormalize norm failure"

  !------------------------------------
  ! spatial product: the m = 0 channel of |psi|^2 with weights sums to
  ! <psi|psi> / sqrt(2 pi)
  !------------------------------------

  allocate (src(Grid_nPoints))
  call Grid_Prolate_SpatialProduct(src, orb1s, orb1s, Grid_Prolate_mmax, Grid_Prolate_mmax, Grid_Prolate_mmax, &
                                   conjgQ_=.true., withWeightsQ_=.true.)

  call check(error, real(sum(src(1:Grid_Prolate_nSpatial)), kind=R64), 1.0_R64 / sqrt(TWOPI), thr=1e-8_R64)
  if (allocated(error)) error stop "T_Grid_Prolate spatial product failure"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "all prolate grid checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
