program T_Utils_SphericalHarmonics
  use M_Utils_Types
  use M_Utils_SphericalHarmonics
  use testdrive, only: check, error_type
  implicit none

  type(error_type), allocatable :: error

  ! Tabulated Gaunt coefficients with an m bound agree with the direct ones,
  ! inside the bound and through the fallback outside of it
  call test_gaunt_table_mbound(error)
  if (allocated(error)) error stop 'SphericalHarmonics Gaunt table test failed'

contains

  subroutine test_gaunt_table_mbound(error)
    type(error_type), allocatable, intent(out) :: error

    integer(I32), parameter :: lmax = 6, mmax = 1
    integer(I32) :: l1, m1, l2, m2, l3
    real(R64) :: maxErr, maxVal

    print *, "Testing Gaunt table with lmax = 6, mmax = 1..."

    call SphericalHarmonics_EnsureGauntTable(lmax, lmax, 2 * lmax, mmax_=mmax)

    call check(error, SphericalHarmonics_gauntTableM .eq. mmax)
    if (allocated(error)) return
    call check(error, size(SphericalHarmonics_gauntTable, 1) .eq. 4 + (lmax - 1) * 3)
    if (allocated(error)) return

    maxErr = 0.0_R64
    maxVal = 0.0_R64
    do l1 = 0, lmax
      do m1 = -l1, l1
        do l2 = 0, lmax
          do m2 = -l2, l2
            do l3 = abs(l1 - l2), l1 + l2
              if (abs(m1 + m2) > l3) cycle
              maxErr = max(maxErr, abs(SphericalHarmonics_GauntTabulated(l1, m1, l2, m2, l3) - &
                                       SphericalHarmonics_GauntCoefficient(l1, m1, l2, m2, l3, m1 + m2)))
              maxVal = max(maxVal, abs(SphericalHarmonics_GauntCoefficient(l1, m1, l2, m2, l3, m1 + m2)))
            end do
          end do
        end do
      end do
    end do

    print *, "Maximum deviation from direct evaluation: ", maxErr, " (max |G| = ", maxVal, ")"

    call check(error, maxErr < 1.0e-14_R64 .and. maxVal > 0.1_R64)

  end subroutine test_gaunt_table_mbound

end program T_Utils_SphericalHarmonics
