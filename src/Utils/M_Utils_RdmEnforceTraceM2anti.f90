! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> For totally antisymmetric tensors (one onefold and one twofold trace).
module M_Utils_RdmEnforceTraceM2anti
  use M_Utils_Types

  implicit none

  ! Expansion coefficients for the unitarily invariant decomposition:
  real(R64), allocatable :: pa
  real(R64), allocatable :: pb

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmEnforceTraceM2anti_Apply(m2anti, onefoldTraces_)

    complex(R64), intent(inout), contiguous            :: m2anti(:, :, :, :)
    complex(R64), intent(inout), contiguous, optional, target :: onefoldTraces_(:, :)

    complex(R64), pointer, contiguous :: onefoldTraces(:, :)
    complex(R64), allocatable, target :: onefoldTracesTmp(:, :)

    integer(I32) :: i1, j1, k
    integer(I32) :: nDim

    nDim = size(m2anti, 1)

    if (present(onefoldTraces_)) then
      onefoldTraces => onefoldTraces_
    else
      allocate (onefoldTracesTmp(nDim, nDim))
      onefoldTracesTmp = 0.0_R64
      onefoldTraces => onefoldTracesTmp
    end if

    do j1 = 1, nDim
      do i1 = 1, nDim
        do k = 1, nDim
          ! Trace: Contract 1st Upper (1) with 1st Lower (3)
          ! For totally antisymmetric tensor, this is the only unique trace.
          onefoldTraces(i1, j1) = onefoldTraces(i1, j1) - m2anti(k, i1, k, j1)
        end do
      end do
    end do

    call AddKroneckerTermsm2anti(m2anti, onefoldTraces)

    if (allocated(onefoldTracesTmp)) deallocate (onefoldTracesTmp)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AddKroneckerTermsm2anti(m2anti, onefoldTraces)

    complex(R64), intent(inout), contiguous :: m2anti(:, :, :, :)
    complex(R64), intent(in), contiguous    :: onefoldTraces(:, :)

    complex(R64) :: twofoldTraces

    integer(I32) :: i1, j1, i2, j2, k
    integer(I32) :: nDim

    nDim = size(m2anti, 1)

    call Setup(nDim)

    twofoldTraces = 0.0_R64
    do k = 1, nDim
      twofoldTraces = twofoldTraces + onefoldTraces(k, k)
    end do

    do j2 = 1, nDim
      do j1 = 1, nDim
        do i2 = 1, nDim
          do i1 = 1, nDim

            ! Twofold trace terms (pa)
            ! delta(i1, j1) * delta(i2, j2) - delta(i1, j2) * delta(i2, j1)
            if (i1 == j1 .and. i2 == j2) m2anti(i1, i2, j1, j2) = m2anti(i1, i2, j1, j2) + pa * twofoldTraces
            if (i1 == j2 .and. i2 == j1) m2anti(i1, i2, j1, j2) = m2anti(i1, i2, j1, j2) - pa * twofoldTraces

            ! Onefold trace terms (pb)
            ! delta(i1, j1)
            if (i1 == j1) m2anti(i1, i2, j1, j2) = m2anti(i1, i2, j1, j2) + pb * onefoldTraces(i2, j2)
            ! -delta(i1, j2)
            if (i1 == j2) m2anti(i1, i2, j1, j2) = m2anti(i1, i2, j1, j2) - pb * onefoldTraces(i2, j1)
            ! -delta(i2, j1)
            if (i2 == j1) m2anti(i1, i2, j1, j2) = m2anti(i1, i2, j1, j2) - pb * onefoldTraces(i1, j2)
            ! delta(i2, j2)
            if (i2 == j2) m2anti(i1, i2, j1, j2) = m2anti(i1, i2, j1, j2) + pb * onefoldTraces(i1, j1)

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup(nDim)

    integer(I32), intent(in) :: nDim
    real(R64) :: x

    if (.not. allocated(pa)) then
      allocate (pa)
      allocate (pb)
    else
      return
    end if

    if (nDim <= 2) Write (*, *) 'WARNING! WARNING! WARNING! WARNING! WARNING! '
    if (nDim <= 2) Write (*, *) 'Contraction consistent Reconstruction is not possible for M2anti at N <= 2!'

    x = dble(nDim)

    ! Coefficients for totally antisymmetric tensor
    pb = 1.0_R64 / (x - 2.0_R64)
    pa = -pb / (x - 1.0_R64)

  end subroutine

end module
