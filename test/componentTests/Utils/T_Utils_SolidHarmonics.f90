!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Real regular solid harmonics (M_Utils_SolidHarmonics, Racah normalization):
!> - closed forms for l <= 3 at a point
!> - pointwise recursion (FillValues) against the polynomial coefficients
!>   (FillCoefficients) up to l = 8
!> - addition theorem 1/|r1 - r2| = sum_lm S_lm(r1) S_lm(r2)/|r2|^(2l+1)
!> - regular at r = 0 (S_00 = 1, all others 0)
program T_Utils_SolidHarmonics
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_SolidHarmonics
  use testdrive, only: check, error_type

  implicit none

  integer(I32), parameter :: lmaxPoly = 8
  integer(I32), parameter :: lmaxAddition = 18

  integer(I32) :: l, m, a, b, c, iPoint, lm
  real(R64) :: r(3), r1(3), r2(3), x, y, z, rr, poly, errClosed, errPoly, scaleMax, sumAddition, errAddition
  real(R64) :: values((lmaxAddition + 1)**2), values2((lmaxAddition + 1)**2)
  real(R64), allocatable :: coeffs(:, :, :, :)
  real(R64), parameter :: points(3, 4) = reshape([0.3_R64, -0.7_R64, 1.1_R64, -1.4_R64, 0.2_R64, 0.5_R64, &
                                                  0.9_R64, 0.9_R64, -0.4_R64, 2.0_R64, -1.5_R64, -0.8_R64], [3, 4])
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  !------------------------------------
  ! closed forms for l <= 3
  !------------------------------------

  r = [0.3_R64, -0.7_R64, 1.1_R64]
  x = r(1)
  y = r(2)
  z = r(3)
  rr = x * x + y * y + z * z

  call SolidHarmonics_FillValues(values, r, 3)

  errClosed = 0.0_R64
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(0, 0)) - 1.0_R64))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(1, -1)) - y))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(1, 0)) - z))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(1, 1)) - x))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(2, -2)) - sqrt(3.0_R64) * x * y))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(2, -1)) - sqrt(3.0_R64) * y * z))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(2, 0)) - 0.5_R64 * (3.0_R64 * z * z - rr)))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(2, 1)) - sqrt(3.0_R64) * x * z))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(2, 2)) - 0.5_R64 * sqrt(3.0_R64) * (x * x - y * y)))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(3, 0)) - 0.5_R64 * z * (5.0_R64 * z * z - 3.0_R64 * rr)))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(3, 2)) - 0.5_R64 * sqrt(15.0_R64) * (x * x - y * y) * z))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(3, -2)) - sqrt(15.0_R64) * x * y * z))
  errClosed = max(errClosed, abs(values(SolidHarmonics_LmIndex(3, 3)) - &
                                 sqrt(10.0_R64) / 4.0_R64 * x * (x * x - 3.0_R64 * y * y)))

  call check(error, errClosed, 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_Utils_SolidHarmonics closed forms failure"

  !------------------------------------
  ! pointwise recursion against the polynomial coefficients
  !------------------------------------

  call SolidHarmonics_FillCoefficients(coeffs, lmaxPoly)

  errPoly = 0.0_R64
  do iPoint = 1, size(points, 2)
    r = points(:, iPoint)
    call SolidHarmonics_FillValues(values, r, lmaxPoly)
    do l = 0, lmaxPoly
      scaleMax = norm2(r)**l
      do m = -l, l
        lm = SolidHarmonics_LmIndex(l, m)
        poly = 0.0_R64
        do c = 0, l
          do b = 0, l - c
            a = l - b - c
            poly = poly + coeffs(a, b, c, lm) * r(1)**a * r(2)**b * r(3)**c
          end do
        end do
        errPoly = max(errPoly, abs(values(lm) - poly) / scaleMax)
      end do
    end do
  end do

  call check(error, errPoly, 0.0_R64, thr=1.0e-13_R64)
  if (allocated(error)) error stop "T_Utils_SolidHarmonics recursion vs polynomial coefficients failure"

  !------------------------------------
  ! addition theorem (|r1|/|r2| ~ 0.1, truncated at l = 18)
  !------------------------------------

  r1 = [0.05_R64, -0.12_R64, 0.08_R64]
  r2 = [0.9_R64, 1.1_R64, -1.3_R64]

  call SolidHarmonics_FillValues(values, r1, lmaxAddition)
  call SolidHarmonics_FillValues(values2, r2, lmaxAddition)

  sumAddition = 0.0_R64
  do l = 0, lmaxAddition
    do m = -l, l
      lm = SolidHarmonics_LmIndex(l, m)
      sumAddition = sumAddition + values(lm) * values2(lm) / norm2(r2)**(2 * l + 1)
    end do
  end do
  errAddition = abs(sumAddition - 1.0_R64 / norm2(r1 - r2)) * norm2(r1 - r2)

  call check(error, errAddition, 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_Utils_SolidHarmonics addition theorem failure"

  !------------------------------------
  ! regular at the origin
  !------------------------------------

  call SolidHarmonics_FillValues(values, [0.0_R64, 0.0_R64, 0.0_R64], 4)

  call check(error, abs(values(1) - 1.0_R64) <= 0.0_R64 .and. maxval(abs(values(2:25))) <= 0.0_R64, &
             "S_00 = 1 and all other S_lm = 0 at the origin")
  if (allocated(error)) error stop "T_Utils_SolidHarmonics origin failure"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(A, ES12.4)') " closed forms l <= 3, max error                = ", errClosed
  write (*, '(A, ES12.4)') " recursion vs polynomial coefficients, l <= 8  = ", errPoly
  write (*, '(A, ES12.4)') " addition theorem, rel. error                  = ", errAddition
  print *
  print *, "all solid harmonics checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
