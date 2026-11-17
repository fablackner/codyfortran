! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Coeffs_Hubbard_PlusSpinSym) S_Coeffs_Hubbard_PlusSpinSym
  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Coeffs_Hubbard_PlusSpinSym_Fabricate
    use M_Utils_Say
    use M_Method_Mb
    use M_Coeffs
    use M_Coeffs_Hubbard

    call Say_Fabricate("coeffs.hubbard.plusSpinSym")

    if (Method_Mb_nBodies(1) .ne. Method_Mb_nBodies(2)) &
      error stop "coeffs.hubbard.plusSpinSym needs equal numbers of spin-up and spin-down particles"
    if (int(Coeffs_Hubbard_nCoeffsUP, I64) * (Coeffs_Hubbard_nCoeffsUP + 1) / 2 > huge(Coeffs_nCoeffs)) &
      error stop "coeffs.hubbard.plusSpinSym: more than 2^31 - 1 coefficients"

    Coeffs_Hubbard_parity = 1
    Coeffs_nCoeffs = int(int(Coeffs_Hubbard_nCoeffsUP, I64) * (Coeffs_Hubbard_nCoeffsUP + 1) / 2, I32)

  end subroutine

end submodule
