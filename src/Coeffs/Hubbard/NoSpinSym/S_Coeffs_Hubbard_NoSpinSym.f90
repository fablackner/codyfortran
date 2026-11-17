! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Coeffs_Hubbard_NoSpinSym) S_Coeffs_Hubbard_NoSpinSym
  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Coeffs_Hubbard_NoSpinSym_Fabricate
    use M_Utils_Say
    use M_Coeffs
    use M_Coeffs_Hubbard

    call Say_Fabricate("coeffs.hubbard.noSpinSym")

    if (int(Coeffs_Hubbard_nCoeffsUP, I64) * Coeffs_Hubbard_nCoeffsDN > huge(Coeffs_nCoeffs)) &
      error stop "coeffs.hubbard.noSpinSym: more than 2^31 - 1 coefficients"

    Coeffs_Hubbard_parity = 0
    Coeffs_nCoeffs = Coeffs_Hubbard_nCoeffsUP * Coeffs_Hubbard_nCoeffsDN

  end subroutine

end submodule
