!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Component test for the Coulomb interaction sysInteraction.fem.coulomb on a
!> mesh graded toward the origin (box +-10, nLocals 6, hMin 0.2, growth 1.6,
!> hMax 2, multipoleLmax 6).
!>
!> Checks:
!> - flags: conjSymmetricQ, timeIndependentQ, bodyTypeIndependentQ
!> - Gaussian orbital exp(-r^2/(2 sigma^2)): its density exp(-r^2/sigma^2) has
!>   the potential pi^1.5 sigma^3 erf(r/sigma)/r (2 pi sigma^2 at the origin
!>   node) and the Hartree integral sqrt(2) pi^2.5 sigma^5
!> - Coulomb integral J(1s, 1s) = 5/8 for the exact hydrogen 1s orbital
!> - swapped pair: ConjugateInteractionPotential(V[a*, b]) = V[b*, a]
!> - Method's bt1_ = bt2_ = -1 are ignored
!> - strength scaling: refabricated with strength 0.5 the potential of the
!>   1s pair halves (to round-off) and J(1s, 1s) becomes 5/16
program T_SysInteraction_Fem_Coulomb
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use M_SysInteraction
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: sigma = 1.0_R64
  integer(I32), parameter :: nRepeat = 5

  integer(I32) :: i, nodeOrigin, iRepeat
  integer(I64) :: tick0, tick1, tickRate
  real(R64) :: errGaussian, errOrigin, hartree, hartreeExact, coulomb1s, errSwap, errSrcSwap, errBt
  real(R64) :: potExact, tFillReal, tFillComplex, errStrength, coulomb1sHalf
  logical :: foundQ
  real(R64), allocatable :: dist(:)
  complex(R64), allocatable :: orb(:), src(:), srcSwap(:), pot(:), potBt(:), potSwap(:), potConj(:), dOrb(:)
  complex(R64), allocatable :: psi1s(:), a(:), b(:), potFull(:), potHalf(:)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/SysInteraction/Fem/Coulomb/T_SysInteraction_Fem_Coulomb.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysInteraction_Fabricate

  call Grid_Setup
  call SysInteraction_Setup

  allocate (dOrb(Grid_nPoints), potConj(Grid_nPoints))

  !==========================================
  call Say_Section("perform action")
  !==========================================

  !------------------------------------
  ! flags
  !------------------------------------

  call check(error, SysInteraction_conjSymmetricQ .and. SysInteraction_timeIndependentQ .and. &
             SysInteraction_bodyTypeIndependentQ, "the Coulomb interaction must set all three flags")
  if (allocated(error)) error stop "T_SysInteraction_Fem_Coulomb flags failure"

  !------------------------------------
  ! Gaussian density: potential and Hartree integral
  !------------------------------------

  dist = Distances([0.0_R64, 0.0_R64, 0.0_R64])
  nodeOrigin = minloc(dist, dim=1)

  orb = cmplx(exp(-dist**2 / (2.0_R64 * sigma**2)), 0.0_R64, kind=R64)

  call SysInteraction_FillInteractionSrc(src, orb, orb)
  call SysInteraction_FillInteractionPotential(pot, src, 0.0_R64)

  call check(error, size(pot) .eq. Grid_nPoints, "potential size differs from Grid_nPoints")
  if (allocated(error)) error stop "T_SysInteraction_Fem_Coulomb potential size failure"

  errGaussian = 0.0_R64
  do i = 1, Grid_nPoints
    if (.not. (dist(i) < 4.0_R64)) cycle
    if (dist(i) > 0.0_R64) then
      potExact = PI**1.5_R64 * sigma**3 * erf(dist(i) / sigma) / dist(i)
    else
      potExact = 2.0_R64 * PI * sigma**2
    end if
    errGaussian = max(errGaussian, abs(real(pot(i), kind=R64) - potExact))
  end do

  errOrigin = abs(real(pot(nodeOrigin), kind=R64) - 2.0_R64 * PI * sigma**2) / (2.0_R64 * PI * sigma**2)

  call check(error, errGaussian, 0.0_R64, thr=1.0e-5_R64)
  if (allocated(error)) error stop "T_SysInteraction_Fem_Coulomb Gaussian potential failure"
  call check(error, errOrigin, 0.0_R64, thr=1.0e-6_R64)
  if (allocated(error)) error stop "T_SysInteraction_Fem_Coulomb Gaussian potential at the origin failure"

  call SysInteraction_MultiplyWithInteractionPotential(dOrb, pot, orb)
  hartree = real(Grid_InnerProduct(orb, dOrb), kind=R64)
  hartreeExact = sqrt(2.0_R64) * PI**2.5_R64 * sigma**5

  call check(error, hartree, hartreeExact, thr=5.0e-7_R64, rel=.true.)
  if (allocated(error)) error stop "T_SysInteraction_Fem_Coulomb Hartree integral failure"

  ! Method passes bt1_ = bt2_ = -1 for body-type-independent interactions
  call SysInteraction_FillInteractionPotential(potBt, src, 1.5_R64, bt1_=-1, bt2_=-1)
  errBt = maxval(abs(potBt - pot))

  call check(error, errBt <= 0.0_R64, "time and body types must not change the potential")
  if (allocated(error)) error stop "T_SysInteraction_Fem_Coulomb body type arguments failure"

  !------------------------------------
  ! Coulomb integral J(1s, 1s) = 5/8 for exp(-r)/sqrt(pi)
  !------------------------------------

  psi1s = cmplx(exp(-dist) / sqrt(PI), 0.0_R64, kind=R64)

  call SysInteraction_FillInteractionSrc(src, psi1s, psi1s)
  call SysInteraction_FillInteractionPotential(pot, src, 0.0_R64)
  call SysInteraction_MultiplyWithInteractionPotential(dOrb, pot, psi1s)
  coulomb1s = real(Grid_InnerProduct(psi1s, dOrb), kind=R64)

  call check(error, coulomb1s, 0.625_R64, thr=5.0e-7_R64)
  if (allocated(error)) error stop "T_SysInteraction_Fem_Coulomb J(1s, 1s) failure"

  !------------------------------------
  ! strength scaling: strength 0.5 halves the potential
  !------------------------------------

  potFull = pot
  call jsonFile % update("sysInteraction.fem.coulomb.strength", 0.5_R64, foundQ)
  if (.not. foundQ) error stop "T_SysInteraction_Fem_Coulomb: sysInteraction.fem.coulomb.strength not found"
  call SysInteraction_Fabricate

  call SysInteraction_FillInteractionPotential(potHalf, src, 0.0_R64)
  errStrength = maxval(abs(potHalf - 0.5_R64 * potFull)) / maxval(abs(potFull))
  call SysInteraction_MultiplyWithInteractionPotential(dOrb, potHalf, psi1s)
  coulomb1sHalf = real(Grid_InnerProduct(psi1s, dOrb), kind=R64)

  call check(error, errStrength, 0.0_R64, thr=1.0e-15_R64)
  if (allocated(error)) error stop "T_SysInteraction_Fem_Coulomb strength scaling failure"
  call check(error, coulomb1sHalf, 0.5_R64 * coulomb1s, thr=1.0e-15_R64)
  if (allocated(error)) error stop "T_SysInteraction_Fem_Coulomb J(1s, 1s) at strength 0.5 failure"

  ! back to strength 1 for the remaining checks
  call jsonFile % update("sysInteraction.fem.coulomb.strength", 1.0_R64, foundQ)
  call SysInteraction_Fabricate

  !------------------------------------
  ! swapped pair: V[b* a] = conjg(V[a* b])
  !------------------------------------

  a = Gaussian(0.8_R64, [0.3_R64, -0.2_R64, 0.5_R64], 0.3_R64)
  b = psi1s * Gaussian(0.1_R64, [-0.4_R64, 0.1_R64, -0.6_R64], -0.7_R64)

  call SysInteraction_FillInteractionSrc(src, a, b)
  call SysInteraction_FillInteractionSrc(srcSwap, b, a)
  errSrcSwap = maxval(abs(srcSwap - conjg(src)))

  call SysInteraction_FillInteractionPotential(pot, src, 0.0_R64)
  call SysInteraction_FillInteractionPotential(potSwap, srcSwap, 0.0_R64)
  call SysInteraction_ConjugateInteractionPotential(potConj, pot)
  errSwap = maxval(abs(potConj - potSwap)) / maxval(abs(potSwap))

  call check(error, errSrcSwap <= 0.0_R64, "the swapped source must be the conjugate")
  if (allocated(error)) error stop "T_SysInteraction_Fem_Coulomb swapped source failure"
  call check(error, errSwap, 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_SysInteraction_Fem_Coulomb ConjugateInteractionPotential failure"

  !------------------------------------
  ! timings (averages over nRepeat calls)
  !------------------------------------

  call SysInteraction_FillInteractionSrc(src, psi1s, psi1s)
  call system_clock(tick0, tickRate)
  do iRepeat = 1, nRepeat
    call SysInteraction_FillInteractionPotential(pot, src, 0.0_R64)
  end do
  call system_clock(tick1)
  tFillReal = real(tick1 - tick0, R64) / real(tickRate, R64) / nRepeat

  call SysInteraction_FillInteractionSrc(src, a, b)
  call system_clock(tick0)
  do iRepeat = 1, nRepeat
    call SysInteraction_FillInteractionPotential(pot, src, 0.0_R64)
  end do
  call system_clock(tick1)
  tFillComplex = real(tick1 - tick0, R64) / real(tickRate, R64) / nRepeat

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(A, I0)') " Grid_nPoints                                     = ", Grid_nPoints
  write (*, '(A, ES12.4)') " Gaussian density, max |V - V_exact| (r < 4)      = ", errGaussian
  write (*, '(A, ES12.4)') " Gaussian density, V(0) vs 2 pi sigma^2, rel.     = ", errOrigin
  write (*, '(A, F18.12, A, F18.12)') " Hartree integral                                 = ", hartree, &
    "  exact ", hartreeExact
  write (*, '(A, F18.12)') " J(1s, 1s) (exact 0.625)                          = ", coulomb1s
  write (*, '(A, ES12.4)') " strength 0.5: max|V - V(1)/2| / max|V(1)|        = ", errStrength
  write (*, '(A, F18.12)') " strength 0.5: J(1s, 1s)                          = ", coulomb1sHalf
  write (*, '(A, ES12.4)') " swapped pair: max|conjg(V_ab) - V_ba| / max|V_ba| = ", errSwap
  write (*, '(A, F9.4, A)') " time FillInteractionPotential (real source)      = ", tFillReal, " s"
  write (*, '(A, F9.4, A)') " time FillInteractionPotential (complex source)   = ", tFillComplex, " s"
  print *
  print *, "all fem interaction checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Distances of all nodes from r0.
  function Distances(r0) result(values)
    real(R64), intent(in) :: r0(3)
    real(R64), allocatable :: values(:)

    integer(I32) :: jx, jy, jz

    allocate (values(Grid_nPoints))
    do jz = 1, Grid_Fem_nz
      do jy = 1, Grid_Fem_ny
        do jx = 1, Grid_Fem_nx
          values(Grid_Fem_Index(jx, jy, jz)) = norm2([Grid_Fem_axes(1) % points(jx), Grid_Fem_axes(2) % points(jy), &
                                                      Grid_Fem_axes(3) % points(jz)] - r0)
        end do
      end do
    end do

  end function

  !> exp(-alpha |r - r0|^2) (1 + i kappa (x - x0)) on the nodes.
  function Gaussian(alpha, r0, kappa) result(values)
    real(R64), intent(in) :: alpha
    real(R64), intent(in) :: r0(3)
    real(R64), intent(in) :: kappa
    complex(R64), allocatable :: values(:)

    integer(I32) :: jx, jy, jz
    real(R64) :: d(3)

    allocate (values(Grid_nPoints))
    do jz = 1, Grid_Fem_nz
      do jy = 1, Grid_Fem_ny
        do jx = 1, Grid_Fem_nx
          d = [Grid_Fem_axes(1) % points(jx), Grid_Fem_axes(2) % points(jy), Grid_Fem_axes(3) % points(jz)] - r0
          values(Grid_Fem_Index(jx, jy, jz)) = exp(-alpha * sum(d**2)) * cmplx(1.0_R64, kappa * d(1), kind=R64)
        end do
      end do
    end do

  end function

end program
