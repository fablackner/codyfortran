! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Two-scan implementation of radial Poisson solver.
submodule(M_SysInteraction_Ylm_Coulomb_TwoScan) S_SysInteraction_Ylm_Coulomb_TwoScan

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind the two-scan Coulomb solver.
  module subroutine SysInteraction_Ylm_Coulomb_TwoScan_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Ylm

    call Say_Fabricate("sysInteraction.ylm.coulomb.twoScan")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Ylm_FillInteractionPotentialRadial => FillInteractionPotentialRadial

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Compute radial potential via forward/backward prefix sums.
  !>
  !> Algorithm:
  !>   1. Forward pass: C1(i) = C1(i-1) + ρ(i) × rᵢˡ
  !>   2. Backward pass: C2(i) = C2(i+1) + ρ(i) / rᵢˡ⁺¹
  !>   3. Combine: V(i) = factor × (C1(i)/rᵢˡ⁺¹ + rᵢˡ×C2(i) - ρ(i)/rᵢ)
  !>
  !> The subtraction of ρ(i)/rᵢ corrects for double-counting at j=i.
  !>
  !> @param[out] potLm   Radial potential component Vₗₘ(r)
  !> @param[in]  srcLm   Radial density component ρₗₘ(r) (includes weights)
  !> @param[in]  l       Angular momentum quantum number
  !> @param[in]  m       Magnetic quantum number (unused)
  !> @param[in]  time    Physical time (unused)
  !> @param[in]  bt1_    Target body type (unused)
  !> @param[in]  bt2_    Source body type (unused)
  subroutine FillInteractionPotentialRadial(potLm, srcLm, l, m, time, bt1_, bt2_)
    use M_Utils_Constants, only: PI
    use M_Utils_UnusedVariables
    use M_Grid_Ylm
    use M_SysInteraction_Ylm_Coulomb, only: SysInteraction_Ylm_Coulomb_Strength

    complex(R64), intent(out), contiguous :: potLm(:)
    complex(R64), intent(in), contiguous  :: srcLm(:)
    integer(I32), intent(in)  :: l
    integer(I32), intent(in)  :: m
    real(R64), intent(in)     :: time
    integer(I32), intent(in), optional  :: bt1_
    integer(I32), intent(in), optional  :: bt2_

    integer(I32) :: i, nRad
    real(R64)    :: r, factor
    complex(R64), allocatable :: C1(:), C2(:)
    real(R64)    :: strength

    if (.false.) call UnusedVariables_Mark(m, time, bt1_, bt2_)

    nRad = Grid_Ylm_nRadial
    strength = SysInteraction_Ylm_Coulomb_Strength

    allocate (C1(nRad), C2(nRad))

    ! Forward scan: C1(i) = Σⱼ≤ᵢ ρ(j) × rⱼˡ
    C1(1) = srcLm(1) * Grid_Ylm_radialPoints(1)**l
    do i = 2, nRad
      C1(i) = C1(i - 1) + srcLm(i) * Grid_Ylm_radialPoints(i)**l
    end do

    ! Backward scan: C2(i) = Σⱼ≥ᵢ ρ(j) / rⱼˡ⁺¹
    C2(nRad) = srcLm(nRad) / Grid_Ylm_radialPoints(nRad)**(l + 1)
    do i = nRad - 1, 1, -1
      C2(i) = C2(i + 1) + srcLm(i) / Grid_Ylm_radialPoints(i)**(l + 1)
    end do

    factor = strength * (4.0_R64 * PI / (2.0_R64 * real(l, R64) + 1.0_R64))

    do i = 1, nRad
      r = Grid_Ylm_radialPoints(i)
      ! Combine scans, subtract diagonal to avoid double-counting
      potLm(i) = factor * (C1(i) / r**(l + 1) + (r**l) * C2(i) - srcLm(i) / r)
    end do

    deallocate (C1, C2)
  end subroutine

end submodule
