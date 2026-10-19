!> Contraction-free projections of M_Utils_RdmEnforceTrace* for changing dimensions.
!>
!> Apply without traces projects on the contraction-free kernel, so a second
!> application must not change the result. The expansion coefficients depend on
!> the dimension and are cached; calling with nDim = 5, 6, 5 checks that they are
!> recomputed when the dimension changes.
program T_Utils_RdmEnforceTrace
  use M_Utils_Types
  use M_Utils_RdmEnforceTraceM2ab
  use M_Utils_RdmEnforceTraceM2sym
  use M_Utils_RdmEnforceTraceM2anti
  use M_Utils_RdmEnforceTraceM3aab
  use testdrive, only: check, error_type
  implicit none

  type(error_type), allocatable :: error
  integer(I32), parameter :: dims(3) = [5, 6, 5]
  integer(I32) :: iDim, n
  real(R64), parameter :: tol = 1e-12_R64
  complex(R64), allocatable :: m2(:, :, :, :), m2Once(:, :, :, :), m3(:, :, :, :, :, :), m3Once(:, :, :, :, :, :)

  do iDim = 1, size(dims)
    n = dims(iDim)

    ! alpha-beta 2-RDM block
    call FillM2(m2, n, 0)
    call RdmEnforceTraceM2ab_Apply(m2)
    m2Once = m2
    call RdmEnforceTraceM2ab_Apply(m2)
    print '(A,I3,A,E12.4)', " nDim", n, "  M2ab   |P(Px) - Px|", Diff4(m2, m2Once)
    call check(error, Diff4(m2, m2Once) < tol)
    call Report("M2ab projection idempotent")

    ! symmetric 2-body tensor
    call FillM2(m2, n, 1)
    call RdmEnforceTraceM2sym_Apply(m2)
    m2Once = m2
    call RdmEnforceTraceM2sym_Apply(m2)
    print '(A,I3,A,E12.4)', " nDim", n, "  M2sym  |P(Px) - Px|", Diff4(m2, m2Once)
    call check(error, Diff4(m2, m2Once) < tol)
    call Report("M2sym projection idempotent")

    ! antisymmetric 2-body tensor
    call FillM2(m2, n, -1)
    call RdmEnforceTraceM2anti_Apply(m2)
    m2Once = m2
    call RdmEnforceTraceM2anti_Apply(m2)
    print '(A,I3,A,E12.4)', " nDim", n, "  M2anti |P(Px) - Px|", Diff4(m2, m2Once)
    call check(error, Diff4(m2, m2Once) < tol)
    call Report("M2anti projection idempotent")

    ! alpha-alpha-beta 3-RDM block (antisymmetric in the two alpha indices)
    call FillM3(m3, n)
    call RdmEnforceTraceM3aab_Apply(m3)
    m3Once = m3
    call RdmEnforceTraceM3aab_Apply(m3)
    print '(A,I3,A,E12.4)', " nDim", n, "  M3aab  |P(Px) - Px|", &
      sqrt(sum(abs(m3 - m3Once)**2)) / sqrt(sum(abs(m3Once)**2))
    call check(error, sqrt(sum(abs(m3 - m3Once)**2)) / sqrt(sum(abs(m3Once)**2)) < tol)
    call Report("M3aab projection idempotent")
  end do

contains

  subroutine Report(what)
    character(len=*), intent(in) :: what
    if (allocated(error)) then
      print *, "FAILED: ", what, " at nDim", n
      error stop "T_Utils_RdmEnforceTrace failed"
    end if
  end subroutine

  function Diff4(a, b) result(res)
    complex(R64), intent(in) :: a(:, :, :, :), b(:, :, :, :)
    real(R64) :: res
    res = sqrt(sum(abs(a - b)**2)) / sqrt(sum(abs(b)**2))
  end function

  !> deterministic complex tensor; sym = 1 symmetric, -1 antisymmetric under i1<->i2 and j1<->j2, 0 none
  subroutine FillM2(m, n, sym)
    complex(R64), intent(out), allocatable :: m(:, :, :, :)
    integer(I32), intent(in) :: n
    integer(I32), intent(in) :: sym
    complex(R64), allocatable :: x(:, :, :, :)
    integer(I32) :: i1, i2, j1, j2
    allocate (x(n, n, n, n), m(n, n, n, n))
    do j2 = 1, n
      do j1 = 1, n
        do i2 = 1, n
          do i1 = 1, n
            x(i1, i2, j1, j2) = cmplx(cos(1.3_R64 * i1 + 0.7_R64 * i2 - 0.4_R64 * j1 + 2.1_R64 * j2), &
                                      sin(0.5_R64 * i1 - 1.7_R64 * i2 + 0.9_R64 * j1 + 0.3_R64 * j2), R64)
          end do
        end do
      end do
    end do
    if (sym .eq. 0) then
      m = x
      return
    end if
    do j2 = 1, n
      do j1 = 1, n
        do i2 = 1, n
          do i1 = 1, n
            m(i1, i2, j1, j2) = x(i1, i2, j1, j2) + sym * x(i2, i1, j1, j2) + &
                                sym * x(i1, i2, j2, j1) + x(i2, i1, j2, j1)
          end do
        end do
      end do
    end do
  end subroutine

  subroutine FillM3(m, n)
    complex(R64), intent(out), allocatable :: m(:, :, :, :, :, :)
    integer(I32), intent(in) :: n
    complex(R64), allocatable :: x(:, :, :, :, :, :)
    integer(I32) :: i1, i2, i3, j1, j2, j3
    allocate (x(n, n, n, n, n, n), m(n, n, n, n, n, n))
    do j3 = 1, n
      do j2 = 1, n
        do j1 = 1, n
          do i3 = 1, n
            do i2 = 1, n
              do i1 = 1, n
                x(i1, i2, i3, j1, j2, j3) = cmplx(cos(1.3_R64 * i1 + 0.7_R64 * i2 - 0.4_R64 * i3 + 2.1_R64 * j1 &
                                                      - 0.6_R64 * j2 + 0.2_R64 * j3), &
                                                  sin(0.5_R64 * i1 - 1.7_R64 * i2 + 0.9_R64 * i3 + 0.3_R64 * j1 &
                                                      + 1.1_R64 * j2 - 0.8_R64 * j3), R64)
              end do
            end do
          end do
        end do
      end do
    end do
    do j3 = 1, n
      do j2 = 1, n
        do j1 = 1, n
          do i3 = 1, n
            do i2 = 1, n
              do i1 = 1, n
                m(i1, i2, i3, j1, j2, j3) = x(i1, i2, i3, j1, j2, j3) - x(i2, i1, i3, j1, j2, j3) - &
                                            x(i1, i2, i3, j2, j1, j3) + x(i2, i1, i3, j2, j1, j3)
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine

end program T_Utils_RdmEnforceTrace
