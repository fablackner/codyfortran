program T_Utils_Legendre
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Legendre
  use testdrive, only: check, error_type

  implicit none

  integer(I32) :: l, m, tau
  real(R64) :: x, y, wronskian, heineSum
  real(R64) :: P(0:40), Q(0:40), Py(0:40), Qy(0:40)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  !------------------------------------
  ! closed forms of P on the cut
  !------------------------------------

  x = 0.3_R64
  call Legendre_FillP(P, 3, 0, x)
  call check(error, P(2), 0.5_R64 * (3.0_R64 * x * x - 1.0_R64), thr=1e-14_R64)
  if (allocated(error)) error stop "P_2^0 on cut failure"

  x = 0.5_R64
  call Legendre_FillP(P, 2, 1, x)
  call check(error, P(2), -3.0_R64 * x * sqrt(1.0_R64 - x * x), thr=1e-14_R64)
  if (allocated(error)) error stop "P_2^1 on cut failure"

  x = -0.4_R64
  call Legendre_FillP(P, 3, 2, x)
  call check(error, P(3), 15.0_R64 * x * (1.0_R64 - x * x), thr=1e-13_R64)
  if (allocated(error)) error stop "P_3^2 on cut failure"

  !------------------------------------
  ! closed forms of P off the cut
  !------------------------------------

  x = 2.0_R64
  call Legendre_FillP(P, 2, 0, x)
  call check(error, P(2), 0.5_R64 * (3.0_R64 * x * x - 1.0_R64), thr=1e-13_R64)
  if (allocated(error)) error stop "P_2^0 off cut failure"

  x = 3.0_R64
  call Legendre_FillP(P, 2, 2, x)
  call check(error, P(2), 3.0_R64 * (x * x - 1.0_R64), thr=1e-12_R64)
  if (allocated(error)) error stop "P_2^2 off cut failure"

  !------------------------------------
  ! closed forms of Q off the cut
  !------------------------------------

  x = 2.0_R64
  call Legendre_FillQ(Q, 2, 0, x)
  call check(error, Q(0), 0.5_R64 * log((x + 1.0_R64) / (x - 1.0_R64)), thr=1e-14_R64)
  if (allocated(error)) error stop "Q_0 failure"
  call check(error, Q(1), x * Q(0) - 1.0_R64, thr=1e-14_R64)
  if (allocated(error)) error stop "Q_1 failure"
  call check(error, Q(2), 0.5_R64 * (3.0_R64 * x * x - 1.0_R64) * Q(0) - 1.5_R64 * x, thr=1e-13_R64)
  if (allocated(error)) error stop "Q_2 failure"

  call Legendre_FillQ(Q, 1, 1, x)
  call check(error, Q(0), -1.0_R64 / sqrt(x * x - 1.0_R64), thr=1e-14_R64)
  if (allocated(error)) error stop "Q_0^1 failure"
  call check(error, Q(1), &
             sqrt(x * x - 1.0_R64) * 0.5_R64 * log((x + 1.0_R64) / (x - 1.0_R64)) - x / sqrt(x * x - 1.0_R64), &
             thr=1e-13_R64)
  if (allocated(error)) error stop "Q_1^1 failure"

  !------------------------------------
  ! Wronskian identity (validates the Green's-function normalization):
  ! (x^2-1)(P_l^m Q_l^m' - P_l^m' Q_l^m) = (l+m)(P_(l-1)^m Q_l^m - P_l^m Q_(l-1)^m)
  !                                      = (-1)^(m+1) (l+m)!/(l-m)!
  !------------------------------------

  do m = 0, 3
    do l = max(1, m), 8
      do tau = 1, 3
        x = 1.0_R64 + 0.5_R64 * tau  ! x = 1.5, 2.0, 2.5
        call Legendre_FillP(P, l, m, x)
        call Legendre_FillQ(Q, l, m, x)
        wronskian = (l + m) * (P(l - 1) * Q(l) - P(l) * Q(l - 1))
        call check(error, wronskian, (-1.0_R64)**(m + 1) * Legendre_FactorialRatio(l, m), &
                   thr=1e-10_R64 * Legendre_FactorialRatio(l, m) + 1e-12_R64)
        if (allocated(error)) error stop "Wronskian identity failure"
      end do
    end do
  end do

  !------------------------------------
  ! Heine identity: sum_tau (2 tau + 1) P_tau(x) Q_tau(y) = 1/(y-x) for y > x > 1
  !------------------------------------

  x = 1.05_R64
  y = 2.0_R64
  call Legendre_FillP(Py, 40, 0, x)
  call Legendre_FillQ(Qy, 40, 0, y)
  heineSum = 0.0_R64
  do tau = 0, 40
    heineSum = heineSum + (2 * tau + 1) * Py(tau) * Qy(tau)
  end do
  call check(error, heineSum, 1.0_R64 / (y - x), thr=1e-10_R64)
  if (allocated(error)) error stop "Heine identity failure"

  !------------------------------------
  ! normalization factor: NormFactorP(l,m)^2 * (l+m)!/(l-m)! = (2l+1)/2
  !------------------------------------

  call check(error, Legendre_NormFactorP(5, 3)**2 * Legendre_FactorialRatio(5, 3), &
             11.0_R64 / 2.0_R64, thr=1e-12_R64)
  if (allocated(error)) error stop "NormFactorP failure"

  print *
  print *, "all Legendre checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
