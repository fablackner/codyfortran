!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The eigen variant of the exact field-free part (fieldFree.eigen,
!> S_FieldFree_Eigen) on its own: H atom on a Ylm FEDVR grid (lmax = mmax = 2,
!> rmax 20, 10 elements of 8 points), h0 = T - 1/r, three stored orbitals.
!>
!> Checks without a cut, with three complex orbitals that fill every (l, m)
!> channel:
!> - unitarity in the grid metric: the overlap matrix is unchanged by
!>   exp(-i h0 dt)
!> - exp(-i h0 dt) exp(-i h0 dt) = exp(-2i h0 dt) and
!>   exp(+i h0 dt) exp(-i h0 dt) = 1
!> - FieldFree_ApplyPhi of one orbital set (phi_0) is FieldFree_Propagate
!> Checks with the lowest eigenvector u_l of the blocks l = 0, 1, 2 (channels
!> (0, 0), (1, 1), (2, -2)), obtained in imaginary time,
!> u_l ~ exp(-h0 tau) of an orbital in that channel:
!> - h0 u_l = E_l u_l with h0 applied here through SysKinetic and SysPotential
!>   (the blocks of the module are those of the operators)
!> - exp(-i h0 dt) u_l = exp(-i E_l dt) u_l
!> - phi_0 + phi_1 of FieldFree_ApplyPhi on (u_l, u_l):
!>   (exp(z) + (exp(z) - 1)/z) u_l, z = -i E_l dt
!> Checks with a spectral cut (fieldFree.eigen.cut = program.cut, new
!> fabricate and setup): P = FieldFree_ProjectBelowCut is a projector that
!> removes weight, P P = P, the propagation maps into its range and ignores
!> the removed part, exp(-i h0 dt) P = P exp(-i h0 dt) P.
program T_FieldFree_EigenYlm
  use M_Utils_Types
  use M_Utils_Constants
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Ylm
  use M_SysKinetic
  use M_SysPotential
  use M_FieldFree
  use M_Method
  use M_Orbs
  use testdrive, only: check, error_type

  implicit none

  !> Stored orbitals (method.mb.orbBased.nOrbs) and their channels (l, m)
  integer(I32), parameter :: N_ORBS = 3
  integer(I32), parameter :: L_OF_ORB(N_ORBS) = [0, 1, 2]
  integer(I32), parameter :: M_OF_ORB(N_ORBS) = [0, 1, -2]

  type(error_type), allocatable :: error
  complex(R64), allocatable :: phi(:, :), work(:, :), work2(:, :), eigenOrbs(:, :), hOrb(:), list(:, :, :)
  complex(R64), allocatable :: externalPotential(:)
  complex(R64) :: z
  real(R64) :: dt, dtau, cut, energy(N_ORBS), dev, removed
  integer(I32) :: nImagSteps, iOrb, iStep
  logical :: foundQ

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/FieldFree/T_FieldFree_EigenYlm.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call FieldFree_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call FieldFree_Setup

  call Say_Fabricate("program")
  dt = Json_Get("program.timeStep", 0.05_R64)
  dtau = Json_Get("program.imagTimeStep", 200.0_R64)
  nImagSteps = Json_Get("program.nImagTimeSteps", 10)
  cut = Json_Get("program.cut", 50.0_R64)

  if (Orbs_nOrbsInState .ne. N_ORBS) error stop "T_FieldFree_EigenYlm: the JSON must set nOrbs = [3]"
  call check(error, associated(FieldFree_Propagate) .and. associated(FieldFree_ApplyPhi) .and. .not. FieldFree_cutQ)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: fieldFree.eigen did not bind its procedures"

  !==========================================
  call Say_Section("perform action")
  !==========================================

  call SysPotential_FillExternalPotential(externalPotential, 0.0_R64)
  allocate (work(Grid_nPoints, N_ORBS), work2(Grid_nPoints, N_ORBS), hOrb(Grid_nPoints))
  call FillTestOrbitals(phi)

  ! 1. Unitarity
  work(:, :) = phi(:, :)
  call FieldFree_Propagate(work, dt, .false.)
  dev = maxval(abs(OverlapMatrix(work) - OverlapMatrix(phi)))
  print *, "overlap matrix after exp(-i h0 dt), deviation:        ", dev
  call check(error, dev < 1e-12_R64)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: the propagation is not unitary"

  ! 2. Group property and inverse
  work2(:, :) = phi(:, :)
  call FieldFree_Propagate(work2, 2.0_R64 * dt, .false.)
  call FieldFree_Propagate(work, dt, .false.)
  dev = MaxNorm(work - work2)
  print *, "exp(-i h0 dt)^2 - exp(-2i h0 dt):                     ", dev
  call check(error, dev < 1e-11_R64)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: two steps of dt differ from one step of 2 dt"

  call FieldFree_Propagate(work, -2.0_R64 * dt, .false.)
  dev = MaxNorm(work - phi)
  print *, "exp(+2i h0 dt) exp(-2i h0 dt) - 1:                    ", dev
  call check(error, dev < 1e-11_R64)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: the backward step is not the inverse"

  ! 3. phi_0 is the propagation
  allocate (list(Grid_nPoints, N_ORBS, 0:0))
  list(:, :, 0) = phi(:, :)
  call FieldFree_ApplyPhi(work, list, dt)
  work2(:, :) = phi(:, :)
  call FieldFree_Propagate(work2, dt, .false.)
  dev = MaxNorm(work - work2)
  print *, "phi_0(-i h0 dt) - exp(-i h0 dt):                      ", dev
  call check(error, dev < 1e-12_R64)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: phi_0 is not the propagation"
  deallocate (list)

  ! 4. Lowest eigenvector of every l block from imaginary time
  call FillChannelOrbitals(eigenOrbs)
  do iStep = 1, nImagSteps
    call FieldFree_Propagate(eigenOrbs, dtau, .true.)
    do iOrb = 1, N_ORBS
      eigenOrbs(:, iOrb) = eigenOrbs(:, iOrb) / MetricNorm(eigenOrbs(:, iOrb))
    end do
  end do

  dev = 0.0_R64
  do iOrb = 1, N_ORBS
    call ApplyH0(hOrb, eigenOrbs(:, iOrb))
    energy(iOrb) = real(Grid_InnerProduct(eigenOrbs(:, iOrb), hOrb), kind=R64)
    dev = max(dev, MetricNorm(hOrb - energy(iOrb) * eigenOrbs(:, iOrb)))
  end do
  print *, "lowest eigenvalues of the blocks l = 0, 1, 2:         ", energy
  print *, "|h0 u - E u| with the operators of SysKinetic and SysPotential: ", dev
  call check(error, dev < 1e-8_R64)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: exp(-h0 tau) does not converge to an eigenvector of h0"

  ! The 1s and 2p levels of hydrogen fit into the box; the 3d level (-1/18)
  ! is lowered by the free boundary at rmax = 20 (-0.0571, A_Grid.md)
  call check(error, abs(energy(1) + 0.5_R64) < 1e-6_R64 .and. abs(energy(2) + 0.125_R64) < 1e-4_R64 .and. &
             energy(3) > energy(2) .and. energy(3) < 0.0_R64)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: the hydrogen levels are wrong"

  work(:, :) = eigenOrbs(:, :)
  call FieldFree_Propagate(work, dt, .false.)
  dev = 0.0_R64
  do iOrb = 1, N_ORBS
    dev = max(dev, MetricNorm(work(:, iOrb) - exp(-IU * energy(iOrb) * dt) * eigenOrbs(:, iOrb)))
  end do
  print *, "exp(-i h0 dt) u - exp(-i E dt) u:                     ", dev
  call check(error, dev < 1e-10_R64)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: wrong phase of an eigenvector"

  allocate (list(Grid_nPoints, N_ORBS, 0:1))
  list(:, :, 0) = eigenOrbs(:, :)
  list(:, :, 1) = eigenOrbs(:, :)
  call FieldFree_ApplyPhi(work, list, dt)
  dev = 0.0_R64
  do iOrb = 1, N_ORBS
    z = -IU * energy(iOrb) * dt
    dev = max(dev, MetricNorm(work(:, iOrb) - (exp(z) + (exp(z) - 1.0_R64) / z) * eigenOrbs(:, iOrb)))
  end do
  print *, "(phi_0 + phi_1) u - (exp(z) + (exp(z) - 1)/z) u:      ", dev
  call check(error, dev < 1e-10_R64)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: wrong phi_1 of an eigenvector"
  deallocate (list)

  ! 5. Spectral cut: read by a new fabricate
  call jsonFile % update("fieldFree.eigen.cut", cut, foundQ)
  if (.not. foundQ) error stop "T_FieldFree_EigenYlm: fieldFree.eigen.cut missing in the JSON"
  call FieldFree_Fabricate
  call FieldFree_Setup
  call check(error, FieldFree_cutQ)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: FieldFree_cutQ is not set with a cut"

  work(:, :) = phi(:, :)
  call FieldFree_ProjectBelowCut(work)
  removed = MaxNorm(work - phi)
  work2(:, :) = work(:, :)
  call FieldFree_ProjectBelowCut(work2)
  dev = MaxNorm(work2 - work)
  print *, "cut", cut, ": |P phi - phi| =", removed, ", |P P phi - P phi| =", dev
  ! The smooth test orbitals have little weight above the cut (9e-5 at cut 50)
  call check(error, removed > 1e-6_R64 .and. dev < 1e-12_R64)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: FieldFree_ProjectBelowCut is not a projector that removes weight"

  ! exp(-i h0 dt) phi = exp(-i h0 dt) P phi, and the result is in the range of P
  work2(:, :) = phi(:, :)
  call FieldFree_Propagate(work2, dt, .false.)
  call FieldFree_Propagate(work, dt, .false.)
  dev = MaxNorm(work2 - work)
  work(:, :) = work2(:, :)
  call FieldFree_ProjectBelowCut(work)
  print *, "propagation with the cut: |U phi - U P phi| =", dev, ", |P U phi - U phi| =", MaxNorm(work - work2)
  call check(error, dev < 1e-12_R64 .and. MaxNorm(work - work2) < 1e-12_R64)
  if (allocated(error)) error stop "T_FieldFree_EigenYlm failure: the propagation with a cut leaves the kept space"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Three complex orbitals with weight in every (l, m) channel, normalized in
  !> the grid metric: r^l exp(-a r) with a channel-dependent complex factor.
  subroutine FillTestOrbitals(orbs)
    complex(R64), intent(out), allocatable :: orbs(:, :)

    real(R64) :: r, l, m
    integer(I32) :: iGrid, i

    allocate (orbs(Grid_nPoints, N_ORBS))
    do iGrid = 1, Grid_nPoints
      r = Grid_Ylm_rCoord(iGrid)
      l = real(Grid_Ylm_lCoord(iGrid), kind=R64)
      m = real(Grid_Ylm_mCoord(iGrid), kind=R64)
      orbs(iGrid, 1) = r**nint(l) * exp(-r) * cmplx(1.0_R64 + l / 5.0_R64, m / 3.0_R64, kind=R64)
      orbs(iGrid, 2) = r**nint(l) * exp(-r / 2.0_R64) * cmplx(1.0_R64, -m / 2.0_R64, kind=R64)
      orbs(iGrid, 3) = r**(nint(l) + 1) * exp(-r) * cmplx(l - 1.0_R64, 1.0_R64 + m, kind=R64)
    end do
    do i = 1, N_ORBS
      orbs(:, i) = orbs(:, i) / MetricNorm(orbs(:, i))
    end do

  end subroutine

  !> Orbital i lives in the channel (L_OF_ORB(i), M_OF_ORB(i)) only:
  !> r^l exp(-r/(l + 1)) there.
  subroutine FillChannelOrbitals(orbs)
    complex(R64), intent(out), allocatable :: orbs(:, :)

    real(R64) :: r
    integer(I32) :: iGrid, i

    allocate (orbs(Grid_nPoints, N_ORBS))
    orbs(:, :) = (0.0_R64, 0.0_R64)
    do i = 1, N_ORBS
      do iGrid = 1, Grid_nPoints
        if (Grid_Ylm_lCoord(iGrid) .ne. L_OF_ORB(i) .or. Grid_Ylm_mCoord(iGrid) .ne. M_OF_ORB(i)) cycle
        r = Grid_Ylm_rCoord(iGrid)
        orbs(iGrid, i) = r**L_OF_ORB(i) * exp(-r / real(L_OF_ORB(i) + 1, kind=R64))
      end do
    end do

  end subroutine

  !> hOrb = (T + V_ext) orb with the operators of SysKinetic and SysPotential.
  subroutine ApplyH0(hOrb, orb)
    complex(R64), intent(out), contiguous :: hOrb(:)
    complex(R64), intent(in), contiguous :: orb(:)

    complex(R64), allocatable :: vOrb(:)

    allocate (vOrb(size(orb)))
    call SysKinetic_MultiplyWithKineticOp(hOrb, orb, 0.0_R64, bt_=1)
    call SysPotential_MultiplyWithExternalPotential(vOrb, externalPotential, orb)
    hOrb(:) = hOrb(:) + vOrb(:)

  end subroutine

  !> Norm of an orbital in the grid metric.
  function MetricNorm(orb) result(res)
    complex(R64), intent(in), contiguous :: orb(:)
    real(R64) :: res

    res = sqrt(real(Grid_InnerProduct(orb, orb), kind=R64))

  end function

  !> Largest metric norm of the columns.
  function MaxNorm(orbs) result(res)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    real(R64) :: res

    integer(I32) :: i

    res = 0.0_R64
    do i = 1, size(orbs, 2)
      res = max(res, MetricNorm(orbs(:, i)))
    end do

  end function

  !> Overlap matrix <orbs_i|orbs_j> in the grid metric.
  function OverlapMatrix(orbs) result(res)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64) :: res(size(orbs, 2), size(orbs, 2))

    integer(I32) :: i, j

    do j = 1, size(orbs, 2)
      do i = 1, size(orbs, 2)
        res(i, j) = Grid_InnerProduct(orbs(:, i), orbs(:, j))
      end do
    end do

  end function

end program
