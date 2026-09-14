! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module M_Utils_RdmProjection
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmProjection_Apply(m2ab, n2ab)
    complex(R64), intent(inout), contiguous :: m2ab(:, :, :, :)
    complex(R64), intent(in), contiguous    :: n2ab(:, :, :, :)

    complex(R64) :: overlap
    real(R64)    :: normSq

    ! Hilbert-Schmidt inner product: <A, B> = tr(A^dag B) = sum(conjg(A_ij) * B_ij)
    normSq = real(sum(conjg(n2ab) * n2ab), R64)

    if (normSq > epsilon(1.0_R64)) then
      overlap = sum(conjg(n2ab) * m2ab)
      m2ab = m2ab - (overlap / normSq) * n2ab
    end if

  end subroutine

end module
