! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for the prolate LCAO orbital initializer.
submodule(M_OrbsInit_Prolate_Lcao) S_OrbsInit_Prolate_Lcao

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Read the LCAO quantum-number arrays from JSON and bind InitializeOrb.
  module subroutine OrbsInit_Prolate_Lcao_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit

    call Say_Fabricate("orbsInit.prolate.lcao")

    OrbsInit_Prolate_Lcao_charge = Json_Get("orbsInit.prolate.lcao.charge", 1.0_R64)
    OrbsInit_Prolate_Lcao_n = Json_Get("orbsInit.prolate.lcao.n", [1])
    OrbsInit_Prolate_Lcao_l = Json_Get("orbsInit.prolate.lcao.l", [0])
    OrbsInit_Prolate_Lcao_m = Json_Get("orbsInit.prolate.lcao.m", [0])
    OrbsInit_Prolate_Lcao_geradeQ = Json_Get("orbsInit.prolate.lcao.geradeQ", [.true.])

    OrbsInit_InitializeOrb => InitializeOrb

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Fill one orbital with the gerade/ungerade LCAO combination.
!>
!> @details
!> The atomic orbital R_nl(r) Theta_lm(cos th) is evaluated at both foci and
!> combined with relative sign s = +-(-1)^l so that the result has the
!> requested inversion parity. The value lands in the azimuthal channel m of
!> the orbital; all other channels are zero. For m /= 0 the xi = 1 axis point
!> is zeroed (the channel Dirichlet condition of the prolate grid). The
!> framework orthonormalizes the orbitals afterwards.
!>
!> @param[out] orb  Orbital on the full prolate grid (flattened channel order).
!> @param[in]  ind  Orbital index selecting (n, l, m, geradeQ) from the config.
!> @param[in]  bt_  Body type (unused; the guess is species-independent).
  subroutine InitializeOrb(orb, ind, bt_)
    use M_Utils_UnusedVariables
    use M_Grid_Prolate

    complex(R64), intent(out), contiguous :: orb(:)
    integer(I32), intent(in)              :: ind
    integer(I32), intent(in), optional    :: bt_

    complex(R64), allocatable :: fM(:)
    real(R64) :: xi, eta, signFactor, r1, r2, cosTheta1, cosTheta2
    integer(I32) :: n, l, m, iXi, iEta, iSpatial

    if (.false.) call UnusedVariables_Mark(bt_)

    n = OrbsInit_Prolate_Lcao_n(ind)
    l = OrbsInit_Prolate_Lcao_l(ind)
    m = OrbsInit_Prolate_Lcao_m(ind)

    if (abs(m) > Grid_Prolate_mmax) error stop "orbsInit.prolate.lcao: |m| exceeds grid mmax"
    if (l < abs(m) .or. n - 1 < l) error stop "orbsInit.prolate.lcao: need |m| <= l <= n-1"

    if (OrbsInit_Prolate_Lcao_geradeQ(ind)) then
      signFactor = (-1.0_R64)**l
    else
      signFactor = -(-1.0_R64)**l
    end if

    allocate (fM(Grid_Prolate_nSpatial))

    do iEta = 1, Grid_Prolate_nEta
      eta = Grid_Prolate_etaPoints(iEta)
      do iXi = 1, Grid_Prolate_nXi
        xi = Grid_Prolate_xiPoints(iXi)
        iSpatial = (iEta - 1) * Grid_Prolate_nXi + iXi

        ! distances and polar angles relative to the foci at z = -a (r1)
        ! and z = +a (r2); eta nodes are interior, so both denominators are
        ! bounded away from zero
        r1 = Grid_Prolate_a * (xi + eta)
        cosTheta1 = (1.0_R64 + xi * eta) / (xi + eta)
        r2 = Grid_Prolate_a * (xi - eta)
        cosTheta2 = (xi * eta - 1.0_R64) / (xi - eta)

        fM(iSpatial) = RadialFunction(r1, n, l) * ThetaFunction(cosTheta1, l, abs(m)) &
                       + signFactor * RadialFunction(r2, n, l) * ThetaFunction(cosTheta2, l, abs(m))

        ! channel Dirichlet condition on the internuclear axis
        if (m .ne. 0 .and. iXi .eq. 1) fM(iSpatial) = (0.0_R64, 0.0_R64)

      end do
    end do

    orb = (0.0_R64, 0.0_R64)
    call Grid_Prolate_AddMComponent(orb, m, fM)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Hydrogen-like radial function R_nl(r) up to normalization.
!>
!> @details
!> R_nl(r) = rho^l L_(n-l-1)^(2l+1)(rho) exp(-rho/2) with rho = 2 Z r / n.
  function RadialFunction(r, n, l) result(res)
    use M_Utils_SfGslLib

    real(R64)                :: res
    real(R64), intent(in)    :: r
    integer(I32), intent(in) :: n, l

    real(R64) :: rho

    rho = 2.0_R64 * OrbsInit_Prolate_Lcao_charge * r / n
    res = exp(-rho / 2.0_R64) * rho**l * SfGslLib_Laguerre(n - l - 1, 2 * l + 1, rho)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Normalized polar factor Theta_lm(x) of the spherical harmonic.
  function ThetaFunction(x, l, mAbs) result(res)
    use M_Utils_Legendre

    real(R64)                :: res
    real(R64), intent(in)    :: x
    integer(I32), intent(in) :: l, mAbs

    real(R64) :: p(0:l)

    call Legendre_FillP(p, l, mAbs, x)
    res = Legendre_NormFactorP(l, mAbs) * p(l)

  end function

end submodule
