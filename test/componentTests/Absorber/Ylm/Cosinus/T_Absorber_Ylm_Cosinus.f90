!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Verifies the radial cosine absorber on the Ylm grid.
!>
!> Applies the absorber to constant orbitals and checks that the imprinted
!> mask matches the analytic form M(r) = cos^(1/order)(π/2 · ξ(r)) with
!> ξ(r) = (r - onset)/(r_max - onset): unity for r < onset, zero at the
!> outermost radial point, and identical for every (l,m) channel.
program T_Absorber_Ylm_Cosinus
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Ylm
  use M_Orbs
  use M_Absorber
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: iGrid, order
  real(R64)    :: r, onset, rLast, t, maskExpected, maxError
  complex(R64), allocatable :: orbs(:, :)
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/Absorber/Ylm/Cosinus/T_Absorber_Ylm_Cosinus.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Absorber_Fabricate

  call Grid_Setup
  call Absorber_Setup

  onset = Json_Get("absorber.ylm.cosinus.onset", 0.0_R64)
  order = Json_Get("absorber.ylm.cosinus.order", 6)
  rLast = Grid_Ylm_radialPoints(Grid_Ylm_nRadial)

  ! Constant orbitals: after applying the absorber they equal the mask itself.
  Orbs_nOrbsInState = 2
  allocate (orbs(Grid_nPoints, Orbs_nOrbsInState))
  orbs(:, 1) = (1.0_R64, 0.0_R64)
  orbs(:, 2) = (0.0_R64, 2.0_R64)

  call Absorber_ApplyAbsorber(orbs)

  !==========================================
  call Say_Section("check mask against analytic form")
  !==========================================

  maxError = 0.0_R64
  do iGrid = 1, Grid_nPoints
    r = Grid_Ylm_rCoord(iGrid)

    if (r < onset) then
      maskExpected = 1.0_R64
    else
      t = (r - onset) / (rLast - onset)
      maskExpected = cos(0.5_R64 * PI * t)**(1.0_R64 / order)
    end if

    maxError = max(maxError, abs(real(orbs(iGrid, 1)) - maskExpected))
    maxError = max(maxError, abs(aimag(orbs(iGrid, 2)) - 2.0_R64 * maskExpected))
  end do

  write (*, '(A, E12.5)') "max deviation from analytic mask: ", maxError

  ! Mask (applied to both orbitals) must match the analytic form on every
  ! (r,l,m) point — this also proves all (l,m) channels see the same M(r).
  call check(error, maxError < 1e-14_R64)
  if (allocated(error)) error stop "T_Absorber_Ylm_Cosinus failure: mask off analytic form"

  ! The mask must vanish at the outermost radial point of every channel.
  ! Floating-point cos(π/2) ~ 1e-16 is lifted to ~1e-4 by the 1/order
  ! exponent, so "vanish" means small, not zero.
  call check(error, maxval(abs(orbs(Grid_Ylm_nRadial::Grid_Ylm_nRadial, 1))) < 1e-3_R64)
  if (allocated(error)) error stop "T_Absorber_Ylm_Cosinus failure: mask nonzero at outer boundary"

  deallocate (orbs)

  !==========================================
  call Say_Goodbye
  !==========================================

end program
