program T_Grid_Ylm_Nodes
  !-----------------------------------------------------------------------------
  ! Parity-split node transforms (Grid_Ylm_ChannelsToNodes/NodesToChannels)
  ! against the unsplit sums over all Gauss-Legendre nodes with Theta_lm from
  ! Legendre_FillP, for every m (also m < 0) and every channel count nCh (odd
  ! and even) the tables cover; the round trip channels -> nodes -> channels is
  ! the identity (the quadrature is exact up to degree 2 nNodes - 1).
  !-----------------------------------------------------------------------------
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Utils_Legendre
  use M_Grid
  use M_Grid_Ylm
  use stdlib_quadrature, only: gauss_legendre
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: nRad, nNodes, m, nCh, l, q
  real(R64) :: errTo, errFrom, errRound
  real(R64), allocatable :: x(:), w(:), theta(:, :), p(:), re(:, :), im(:, :)
  complex(R64), allocatable :: fM(:, :), fM2(:, :), fNodes(:, :), gNodes(:, :), ref(:, :)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Grid/Ylm/T_Grid_Ylm_Nodes.json", &
                         relativeToProjectDirQ_=.true.)
  call Grid_Fabricate
  call Grid_Setup

  nRad = Grid_Ylm_nRadial
  nNodes = Grid_Ylm_nNodes
  allocate (x(nNodes), w(nNodes), p(0:Grid_Ylm_lmaxNodes))
  call gauss_legendre(x, w)

  !==========================================
  call Say_Section("perform action")
  !==========================================

  errTo = 0.0_R64
  errFrom = 0.0_R64
  errRound = 0.0_R64
  do m = -Grid_Ylm_mmaxNodes, Grid_Ylm_mmaxNodes
    ! Unsplit table theta(l, q) = Theta_lm(x_q), Theta_l,-m = (-1)^m Theta_l|m|
    allocate (theta(abs(m):Grid_Ylm_lmaxNodes, nNodes))
    do q = 1, nNodes
      call Legendre_FillP(p, Grid_Ylm_lmaxNodes, abs(m), x(q))
      do l = abs(m), Grid_Ylm_lmaxNodes
        theta(l, q) = (-1)**(merge(abs(m), 0, m < 0)) * Legendre_NormFactorP(l, abs(m)) * p(l) / sqrt(2.0_R64 * PI)
      end do
    end do

    do nCh = 1, Grid_Ylm_lmaxNodes - abs(m) + 1
      allocate (re(nRad, nCh), im(nRad, nCh), fNodes(nRad, nNodes), fM2(nRad, nCh))
      call random_number(re)
      call random_number(im)
      fM = cmplx(re - 0.5_R64, im - 0.5_R64, kind=R64)
      deallocate (re, im)
      allocate (re(nRad, nNodes), im(nRad, nNodes))
      call random_number(re)
      call random_number(im)
      gNodes = cmplx(re - 0.5_R64, im - 0.5_R64, kind=R64)

      call Grid_Ylm_ChannelsToNodes(fNodes, fM, m)
      ref = matmul(fM, theta(abs(m):abs(m) + nCh - 1, :))
      errTo = max(errTo, maxval(abs(fNodes - ref)))

      call Grid_Ylm_NodesToChannels(fM2, fNodes, m)
      errRound = max(errRound, maxval(abs(fM2 - fM)))

      call Grid_Ylm_NodesToChannels(fM2, gNodes, m)
      ref = matmul(gNodes, transpose(spread(Grid_Ylm_nodeWeights, 1, nCh) * theta(abs(m):abs(m) + nCh - 1, :)))
      errFrom = max(errFrom, maxval(abs(fM2 - ref)))

      deallocate (re, im, fNodes, fM2)
    end do
    deallocate (theta)
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  write (*, '(A, ES12.5)') "max |ChannelsToNodes - unsplit|: ", errTo
  write (*, '(A, ES12.5)') "max |NodesToChannels - unsplit|: ", errFrom
  write (*, '(A, ES12.5)') "max |round trip - identity|:     ", errRound

  call check(error, errTo < 1.0e-12_R64 .and. errFrom < 1.0e-12_R64 .and. errRound < 1.0e-12_R64)
  if (allocated(error)) error stop "T_Grid_Ylm_Nodes failure"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
