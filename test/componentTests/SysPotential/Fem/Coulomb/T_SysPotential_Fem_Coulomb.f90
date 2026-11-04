!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Component test for the static nuclear attraction sysPotential.fem.coulomb.
!>
!> Mesh: box +-8, nLocals 6, growth 1.6, hMax 2, four centers
!>   A (0, 0, 0)     Z = 1, point nucleus (Galerkin), hMin 0.2
!>   B (0, 0, 2.5)   Z = 2, smeared, rc = 0.4,         hMin 0.4
!>   C (1.5, 0, -1)  Z = 3, smeared, rc = 0.3,         hMin 0.4
!>   D (1.5, 0, 0)   Z = 0, ghost (no new mesh lines),  hMin 0.4
!>
!> Checks:
!> - flags (time- and body-type-independent), nuclear repulsion of A, B, C
!> - FillExternalPotential returns the same static potential at any time and
!>   for any body type, and the multiplication is point-wise
!> - erf smearing: static potential = Galerkin(A) + nodal -Z erf(r/rc)/r of
!>   B and C (the ghost adds nothing), including the r = 0 limit
!>   -2Z/(sqrt(pi) rc) on the nucleus nodes of B and C
!> - the point proton A: potential against -1/r away from the nucleus, finite
!>   at the nucleus node (-49.6/hMin on this mesh), and <1s|V_A|1s> = -1 for the exact
!>   hydrogen 1s orbital
program T_SysPotential_Fem_Coulomb
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use M_SysPotential
  use M_SysPotential_Fem_Coulomb
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: centerA(3) = [0.0_R64, 0.0_R64, 0.0_R64]
  real(R64), parameter :: centerB(3) = [0.0_R64, 0.0_R64, 2.5_R64]
  real(R64), parameter :: centerC(3) = [1.5_R64, 0.0_R64, -1.0_R64]
  real(R64), parameter :: chargeB = 2.0_R64, smearingB = 0.4_R64
  real(R64), parameter :: chargeC = 3.0_R64, smearingC = 0.3_R64
  real(R64), parameter :: hMinA = 0.2_R64

  integer(I32) :: i, nodeA, nodeB, nodeC
  real(R64) :: repulsionExact, potScale, errSmeared, errLimitB, errLimitC, potNucleusA
  real(R64) :: limitB, errFunction0, errFunctionSmall, errFunction07
  real(R64) :: errPoint15, errPoint3, energy1s
  real(R64), allocatable :: distA(:), distB(:), distC(:), smeared(:)
  complex(R64), allocatable :: pot0(:), pot1(:), potA(:), potPointA(:), orb(:), dOrb(:), psi1s(:)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/SysPotential/Fem/Coulomb/T_SysPotential_Fem_Coulomb.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysPotential_Fabricate

  call Grid_Setup
  call SysPotential_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  !------------------------------------
  ! flags and nuclear repulsion (the ghost D does not enter)
  !------------------------------------

  call check(error, SysPotential_timeIndependentQ .and. SysPotential_bodyTypeIndependentQ, &
             "the nuclear potential must be time- and body-type-independent")
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb flags failure"

  repulsionExact = 1.0_R64 * chargeB / norm2(centerA - centerB) + 1.0_R64 * chargeC / norm2(centerA - centerC) + &
                   chargeB * chargeC / norm2(centerB - centerC)

  call check(error, SysPotential_Fem_Coulomb_nuclearRepulsion, repulsionExact, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb nuclear repulsion failure"

  !------------------------------------
  ! Fill returns the same static potential at any time and for any body type;
  ! pot0 is the static potential checked below
  !------------------------------------

  call SysPotential_FillExternalPotential(pot0, 0.0_R64)
  call SysPotential_FillExternalPotential(pot1, 3.7_R64, bt_=2)

  call check(error, size(pot0) .eq. Grid_nPoints .and. size(pot1) .eq. Grid_nPoints, "potential size differs from Grid_nPoints")
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb Fill size failure"
  call check(error, maxval(abs(pot1 - pot0)) <= 0.0_R64, "Fill must return the same static potential at any time and body type")
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb Fill copy failure"

  ! point-wise multiplication
  orb = Gaussian(0.7_R64, [0.2_R64, -0.1_R64, 0.4_R64], 0.5_R64)
  allocate (dOrb(Grid_nPoints))
  call SysPotential_MultiplyWithExternalPotential(dOrb, pot0, orb)

  call check(error, maxval(abs(dOrb - pot0 * orb)) <= 0.0_R64, "the multiplication must be point-wise")
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb multiply failure"

  !------------------------------------
  ! smeared nuclei: static = Galerkin(A) + nodal erf potentials of B and C
  !------------------------------------

  distA = Distances(centerA)
  distB = Distances(centerB)
  distC = Distances(centerC)

  allocate (smeared(Grid_nPoints), potA(Grid_nPoints))
  do i = 1, Grid_nPoints
    smeared(i) = ErfPotential(chargeB, smearingB, distB(i)) + ErfPotential(chargeC, smearingC, distC(i))
  end do

  call Grid_Fem_SolvePoissonPointCharges(potA, [-1.0_R64], reshape(centerA, [3, 1]))

  potScale = maxval(abs(pot0))
  errSmeared = maxval(abs(pot0 - potA - smeared)) / potScale

  call check(error, errSmeared, 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb smeared nuclei failure"

  ! the r = 0 limit on the nucleus nodes of B and C
  nodeA = NodeAt(centerA)
  nodeB = NodeAt(centerB)
  nodeC = NodeAt(centerC)

  limitB = -2.0_R64 * chargeB / (sqrt(PI) * smearingB)
  errLimitB = abs(real(pot0(nodeB), kind=R64) - real(potA(nodeB), kind=R64) - &
                  ErfPotential(chargeC, smearingC, distC(nodeB)) - limitB) / abs(limitB)
  errLimitC = abs(real(pot0(nodeC), kind=R64) - real(potA(nodeC), kind=R64) - &
                  ErfPotential(chargeB, smearingB, distB(nodeC)) + 2.0_R64 * chargeC / (sqrt(PI) * smearingC)) / &
              (2.0_R64 * chargeC / (sqrt(PI) * smearingC))

  call check(error, errLimitB, 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb smeared limit (B) failure"
  call check(error, errLimitC, 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb smeared limit (C) failure"

  ! the exported smeared potential: limit, continuity at small r, a value
  errFunction0 = abs(SysPotential_Fem_Coulomb_SmearedPotential(chargeB, smearingB, 0.0_R64) - limitB) / abs(limitB)
  errFunctionSmall = abs(SysPotential_Fem_Coulomb_SmearedPotential(chargeB, smearingB, 1.0e-7_R64) - limitB) / abs(limitB)
  errFunction07 = abs(SysPotential_Fem_Coulomb_SmearedPotential(chargeB, smearingB, 0.7_R64) + &
                      chargeB * erf(0.7_R64 / smearingB) / 0.7_R64)

  call check(error, errFunction0, 0.0_R64, thr=1.0e-15_R64)
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb smeared potential at r = 0 failure"
  call check(error, errFunctionSmall, 0.0_R64, thr=1.0e-12_R64)
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb smeared potential continuity failure"
  call check(error, errFunction07, 0.0_R64, thr=1.0e-15_R64)
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb smeared potential value failure"

  !------------------------------------
  ! the Galerkin point proton A (static minus the smeared nuclei)
  !------------------------------------

  potPointA = pot0 - smeared

  errPoint15 = maxval(abs(real(potPointA, kind=R64) + 1.0_R64 / max(distA, 1.0e-12_R64)), mask=distA > 1.5_R64)
  errPoint3 = maxval(abs(real(potPointA, kind=R64) + 1.0_R64 / max(distA, 1.0e-12_R64)), mask=distA > 3.0_R64)
  potNucleusA = real(potPointA(nodeA), kind=R64)

  call check(error, errPoint15, 0.0_R64, thr=2.0e-5_R64)
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb point proton (r > 1.5) failure"
  call check(error, errPoint3, 0.0_R64, thr=4.0e-6_R64)
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb point proton (r > 3) failure"
  call check(error, potNucleusA < -30.0_R64 / hMinA .and. potNucleusA > -60.0_R64 / hMinA, &
             "the point nucleus must give a finite node value between -60/hMin and -30/hMin")
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb value at the nucleus failure"

  ! <1s|V_A|1s> for the exact hydrogen 1s orbital exp(-r)/sqrt(pi)
  psi1s = cmplx(exp(-distA) / sqrt(PI), 0.0_R64, kind=R64)
  energy1s = real(Grid_InnerProduct(psi1s, potPointA * psi1s), kind=R64)

  call check(error, energy1s, -1.0_R64, thr=5.0e-6_R64)
  if (allocated(error)) error stop "T_SysPotential_Fem_Coulomb <1s|V|1s> failure"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(A, I0)') " Grid_nPoints                                     = ", Grid_nPoints
  write (*, '(A, F20.14)') " nuclear repulsion (exact 4.0397777...)           = ", SysPotential_Fem_Coulomb_nuclearRepulsion
  write (*, '(A, ES12.4)') " static - Galerkin(A) - erf(B, C), rel. max      = ", errSmeared
  write (*, '(A, 2ES12.4)') " r = 0 limit on the nodes of B and C, rel. error = ", errLimitB, errLimitC
  write (*, '(A, 3ES12.4)') " SmearedPotential: r = 0, r = 1e-7, r = 0.7      = ", errFunction0, errFunctionSmall, &
    errFunction07
  write (*, '(A, ES12.4)') " point proton, max |V + 1/r| (r > 1.5)           = ", errPoint15
  write (*, '(A, ES12.4)') " point proton, max |V + 1/r| (r > 3)             = ", errPoint3
  write (*, '(A, F16.10, A, F8.3, A)') " point proton, V at the nucleus node             = ", potNucleusA, &
    "  (", potNucleusA * hMinA, " / hMin)"
  write (*, '(A, F16.12)') " <1s|V_A|1s> (exact -1)                          = ", energy1s
  print *
  print *, "all fem nuclear potential checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> -Z erf(r/rc)/r with the limit -2Z/(sqrt(pi) rc) at r = 0.
  pure function ErfPotential(charge, smearing, r) result(res)
    real(R64), intent(in) :: charge
    real(R64), intent(in) :: smearing
    real(R64), intent(in) :: r
    real(R64) :: res

    if (r > 0.0_R64) then
      res = -charge * erf(r / smearing) / r
    else
      res = -2.0_R64 * charge / (sqrt(PI) * smearing)
    end if

  end function

  !> Distances of all nodes from r0.
  function Distances(r0) result(dist)
    real(R64), intent(in) :: r0(3)
    real(R64), allocatable :: dist(:)

    integer(I32) :: jx, jy, jz

    allocate (dist(Grid_nPoints))
    do jz = 1, Grid_Fem_nz
      do jy = 1, Grid_Fem_ny
        do jx = 1, Grid_Fem_nx
          dist(Grid_Fem_Index(jx, jy, jz)) = norm2([Grid_Fem_axes(1) % points(jx), Grid_Fem_axes(2) % points(jy), &
                                                    Grid_Fem_axes(3) % points(jz)] - r0)
        end do
      end do
    end do

  end function

  !> Index of the node at r0 (stops if r0 is not a node).
  function NodeAt(r0) result(node)
    real(R64), intent(in) :: r0(3)
    integer(I32) :: node

    integer(I32) :: idx(3), iAxis

    do iAxis = 1, 3
      idx(iAxis) = minloc(abs(Grid_Fem_axes(iAxis) % points - r0(iAxis)), dim=1)
      if (abs(Grid_Fem_axes(iAxis) % points(idx(iAxis)) - r0(iAxis)) > 0.0_R64) error stop "NodeAt: not a node"
    end do
    node = Grid_Fem_Index(idx(1), idx(2), idx(3))

  end function

  !> exp(-alpha |r - r0|^2) (1 + i kappa (x - x0)) on the nodes.
  function Gaussian(alpha, r0, kappa) result(field)
    real(R64), intent(in) :: alpha
    real(R64), intent(in) :: r0(3)
    real(R64), intent(in) :: kappa
    complex(R64), allocatable :: field(:)

    integer(I32) :: jx, jy, jz
    real(R64) :: d(3)

    allocate (field(Grid_nPoints))
    do jz = 1, Grid_Fem_nz
      do jy = 1, Grid_Fem_ny
        do jx = 1, Grid_Fem_nx
          d = [Grid_Fem_axes(1) % points(jx), Grid_Fem_axes(2) % points(jy), Grid_Fem_axes(3) % points(jz)] - r0
          field(Grid_Fem_Index(jx, jy, jz)) = exp(-alpha * sum(d**2)) * cmplx(1.0_R64, kappa * d(1), kind=R64)
        end do
      end do
    end do

  end function

end program
