! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> unitarily invariant decomposition utilities for second-order tensors (M2) in RDM theory.
!>
!> Supplies routines to enforce trace conditions and assemble the
!> M2 decomposition with precomputed expansion coefficients.
module M_Utils_RdmEnforceTraceM2ab
  use M_Utils_Types

  implicit none

  ! Expansion coefficients for the unitarily invariant decomposition:
  real(R64), allocatable :: pa(:, :)
  real(R64), allocatable :: pb(:, :)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmEnforceTraceM2ab_Apply(m2ab, onefoldTraces_)

    complex(R64), intent(inout), contiguous            :: m2ab(:, :, :, :)
    complex(R64), intent(inout), contiguous, optional, target :: onefoldTraces_(:, :, :)

    complex(R64), pointer, contiguous :: onefoldTraces(:, :, :)
    complex(R64), allocatable, target :: onefoldTracesTmp(:, :, :)

    integer(I32) :: i1, j1, k
    integer(I32) :: nDim

    nDim = size(m2ab, 1)

    if (present(onefoldTraces_)) then
      onefoldTraces => onefoldTraces_
    else
      allocate (onefoldTracesTmp(nDim, nDim, 2))
      onefoldTracesTmp = 0.0_R64
      onefoldTraces => onefoldTracesTmp
    end if

    do j1 = 1, nDim
      do i1 = 1, nDim
        do k = 1, nDim

          ! Trace 1 (A): Contract 1st Upper (1) with 1st Lower (3)
          onefoldTraces(i1, j1, 1) = onefoldTraces(i1, j1, 1) - m2ab(k, i1, k, j1)

          ! Trace 2 (B): Contract 1st Upper (1) with 2nd Lower (4)
          onefoldTraces(i1, j1, 2) = onefoldTraces(i1, j1, 2) - m2ab(k, i1, j1, k)

        end do
      end do
    end do

    call AddKroneckerTermsM2ab(m2ab, onefoldTraces)

    if (allocated(onefoldTracesTmp)) deallocate (onefoldTracesTmp)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AddKroneckerTermsM2ab(m2ab, onefoldTraces)

    complex(R64), intent(inout), contiguous :: m2ab(:, :, :, :)
    complex(R64), intent(in), contiguous    :: onefoldTraces(:, :, :)

    complex(R64), allocatable :: twofoldTraces(:)

    integer(I32) :: i2, j2, i1, j1, k

    integer(I32) :: nDim

    nDim = size(m2ab, 1)

    call Setup(nDim)

    allocate (twofoldTraces(2))

    twofoldTraces(:) = 0.0_R64
    do k = 1, nDim
      twofoldTraces(1) = twofoldTraces(1) + onefoldTraces(k, k, 1)
      twofoldTraces(2) = twofoldTraces(2) + onefoldTraces(k, k, 2)
    end do

    do j2 = 1, nDim
      do j1 = 1, nDim
        do i2 = 1, nDim
          do i1 = 1, nDim

            ! Perm 1: delta(i1, j1) * delta(i2, j2)
            if (i1 == j1 .and. i2 == j2) then
              m2ab(i1, i2, j1, j2) = m2ab(i1, i2, j1, j2) + &
                                     pa(1, 1) * twofoldTraces(1) + &
                                     pa(1, 2) * twofoldTraces(2)
            end if

            ! Perm 2: delta(i1, j2) * delta(i2, j1)
            if (i1 == j2 .and. i2 == j1) then
              m2ab(i1, i2, j1, j2) = m2ab(i1, i2, j1, j2) + &
                                     pa(2, 1) * twofoldTraces(1) + &
                                     pa(2, 2) * twofoldTraces(2)
            end if

            ! Perm 1: delta(i1, j1)
            if (i1 == j1) then
              m2ab(i1, i2, j1, j2) = m2ab(i1, i2, j1, j2) + &
                                     pb(1, 1) * onefoldTraces(i2, j2, 1) + &
                                     pb(1, 2) * onefoldTraces(i2, j2, 2)
            end if

            ! Perm 2: delta(i1, j2)
            if (i1 == j2) then
              m2ab(i1, i2, j1, j2) = m2ab(i1, i2, j1, j2) + &
                                     pb(2, 1) * onefoldTraces(i2, j1, 1) + &
                                     pb(2, 2) * onefoldTraces(i2, j1, 2)
            end if

            ! Perm 3: delta(i2, j1)
            if (i2 == j1) then
              m2ab(i1, i2, j1, j2) = m2ab(i1, i2, j1, j2) + &
                                     pb(3, 1) * onefoldTraces(i1, j2, 1) + &
                                     pb(3, 2) * onefoldTraces(i1, j2, 2)
            end if

            ! Perm 4: delta(i2, j2)
            if (i2 == j2) then
              m2ab(i1, i2, j1, j2) = m2ab(i1, i2, j1, j2) + &
                                     pb(4, 1) * onefoldTraces(i1, j1, 1) + &
                                     pb(4, 2) * onefoldTraces(i1, j1, 2)
            end if

          end do
        end do
      end do
    end do

    if (allocated(twofoldTraces)) deallocate (twofoldTraces)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup(nDim)

    integer(I32), intent(in) :: nDim
    real(R64) :: x

    if (.not. allocated(pa)) then
      allocate (pa(2, 2))
      allocate (pb(4, 2))
    else
      return
    end if

    if (nDim <= 2) Write (*, *) 'WARNING! WARNING! WARNING! WARNING! WARNING! '
    if (nDim <= 2) Write (*, *) 'Contraction consistent Reconstruction is not possible for M2 at N <= 2!'

    x = dble(nDim)

    pa(1, 1) = -(x**2 + 2.0_R64) / ((x**2 - 4.0_R64) * (x**2 - 1.0_R64))
    pa(1, 2) = (3.0_R64 * x) / ((x**2 - 4.0_R64) * (x**2 - 1.0_R64))

    pa(2, 1) = (3.0_R64 * x) / ((x**2 - 4.0_R64) * (x**2 - 1.0_R64))
    pa(2, 2) = -(x**2 + 2.0_R64) / ((x**2 - 4.0_R64) * (x**2 - 1.0_R64))

    pb(1, 1) = x / (x**2 - 4.0_R64)
    pb(1, 2) = -2.0_R64 / (x**2 - 4.0_R64)

    pb(2, 1) = -2.0_R64 / (x**2 - 4.0_R64)
    pb(2, 2) = x / (x**2 - 4.0_R64)

    pb(3, 1) = -2.0_R64 / (x**2 - 4.0_R64)
    pb(3, 2) = x / (x**2 - 4.0_R64)

    pb(4, 1) = x / (x**2 - 4.0_R64)
    pb(4, 2) = -2.0_R64 / (x**2 - 4.0_R64)

  end subroutine

end module
