! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased) S_Method_Mb_OrbBased

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_NoOpProcedures
    use M_Method
    use M_Method_Mb
    use M_Method_Mb_OrbBased_Tdci
    use M_Method_Mb_OrbBased_Tdhx
    use M_Method_Mb_OrbBased_Mctdhx

    integer(I32) :: i, ibt, index, counter

    call Say_Fabricate("method.mb.orbBased")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    allocate (Method_Mb_OrbBased_nOrbs(Method_Mb_nBodyTypes))
    Method_Mb_OrbBased_nOrbs = Json_Get("method.mb.orbBased.nOrbs", [(0, i=1, Method_Mb_nBodyTypes)])

    Method_Mb_OrbBased_nOrbsSum = sum(Method_Mb_OrbBased_nOrbs)

    allocate (Method_Mb_OrbBased_bodyTypeOfOrb(Method_Mb_OrbBased_nOrbsSum))
    allocate (Method_Mb_OrbBased_nOrbsStart(Method_Mb_nBodyTypes))
    allocate (Method_Mb_OrbBased_nOrbsEnd(Method_Mb_nBodyTypes))

    counter = 0
    do ibt = 1, Method_Mb_nBodyTypes
      Method_Mb_OrbBased_nOrbsStart(ibt) = counter + 1
      do index = 1, Method_Mb_OrbBased_nOrbs(ibt)
        counter = counter + 1
        Method_Mb_OrbBased_bodyTypeOfOrb(counter) = ibt
      end do
      Method_Mb_OrbBased_nOrbsEnd(ibt) = counter
    end do

    Method_Mb_OrbBased_regularizationParameter = Json_Get("method.mb.orbBased.regularizationParameter", 1e-10_R64)

    Method_GetEnergy => GetEnergy
    Method_Mb_OrbBased_FillH1 => FillH1
    Method_Mb_OrbBased_FillH2 => FillH2
    Method_Mb_OrbBased_ApplyKineticOp => ApplyKineticOp
    Method_Mb_OrbBased_ApplyCorrelationOp => ApplyCorrelationOp
    Method_Mb_OrbBased_ApplyHartreeFockOp => ApplyHartreeFockOp
    Method_Mb_OrbBased_ApplySingleBodyOp => ApplySingleBodyOp

    if (Json_GetExistence("sysPotential")) then
      Method_Mb_OrbBased_ApplyPotentialOp => ApplyPotentialOp
    else
      Method_Mb_OrbBased_ApplyPotentialOp => NoOpProcedures_ApplyPotentialOp
    end if

    if (Json_GetExistence("sysInteraction")) then
      Method_Mb_OrbBased_ApplyInteractionOp => ApplyInteractionOp
    else
      Method_Mb_OrbBased_ApplyInteractionOp => NoOpProcedures_ApplyInteractionOp
    end if

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("method.mb.orbBased.tdci")) then
      call Method_Mb_OrbBased_Tdci_Fabricate

    else if (Json_GetExistence("method.mb.orbBased.tdhx")) then
      call Method_Mb_OrbBased_Tdhx_Fabricate

    else if (Json_GetExistence("method.mb.orbBased.mctdhx")) then
      call Method_Mb_OrbBased_Mctdhx_Fabricate

    else
      error stop "method.mb.orbBased is missing one of: tdci, tdhx, mctdhx"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function GetEnergy(time) result(res)
    use M_Utils_RdmObservables
    use M_Method
    use M_Orbs

    real(R64)                :: res
    real(R64), intent(in)    :: time

    call Method_Mb_OrbBased_FillRdm1(Method_Mb_OrbBased_rdm1, Method_state)
    call Method_Mb_OrbBased_FillRdm2(Method_Mb_OrbBased_rdm2, Method_state)
    call Method_Mb_OrbBased_FillH1(Method_Mb_OrbBased_h1, Orbs_orbs, time)
    call Method_Mb_OrbBased_FillH2(Method_Mb_OrbBased_h2, Orbs_orbs, time)

    res = RdmObservables_Energy(Method_Mb_OrbBased_rdm1, &
                                Method_Mb_OrbBased_h1, &
                                Method_Mb_OrbBased_rdm2, &
                                Method_Mb_OrbBased_h2)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH1(h1, orbs, time)
    use M_Method
    use M_Orbs
    use M_Grid
    use M_SysKinetic
    use M_SysPotential

    complex(R64), intent(out), allocatable :: h1(:, :)
    complex(R64), intent(in), contiguous                :: orbs(:, :)
    real(R64), intent(in)                  :: time

    integer(I32) :: nOS, nO, nG
    integer(I32) :: i1, j1, i1bt, j1bt
    complex(R64), allocatable :: orbTmps(:, :)

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    if (.not. allocated(h1)) allocate (h1(nO, nO))
    h1(:, :) = 0.0_R64

    allocate (orbTmps(nG, nOS))
    orbTmps(:, :) = 0.0_R64

    call Method_Mb_OrbBased_ApplyKineticOp(orbTmps, orbs, time)
    call Method_Mb_OrbBased_ApplyPotentialOp(orbTmps, orbs, time)

    do j1 = 1, nOS
      do i1 = 1, nOS

        i1bt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)
        j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

        if (i1bt .ne. j1bt) cycle

        h1(i1, j1) = Grid_InnerProduct(orbs(:, i1), orbTmps(:, j1))

      end do
    end do

    if (Orbs_restrictedQ) then
      h1(nOS + 1:, nOS + 1:) = h1(:nOS, :nOS)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH2(h2, orbs, time)
    use M_Method
    use M_SysInteraction
    use M_Orbs
    use M_Grid

    complex(R64), intent(out), allocatable :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous   :: orbs(:, :)
    real(R64), intent(in)                  :: time

    integer(I32) :: nG, nOS, nO
    integer(I32) :: i1, j1, i2, j2
    integer(I32) :: i1bt, j1bt, i2bt, j2bt
    complex(R64), allocatable :: orbTmps(:, :)

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    if (.not. allocated(h2)) allocate (h2(nO, nO, nO, nO))
    h2(:, :, :, :) = 0.0_R64

    allocate (orbTmps, mold=orbs)

    do j2 = 1, nOS
      do i2 = 1, nOS
        i2bt = Method_Mb_OrbBased_bodyTypeOfOrb(i2)
        j2bt = Method_Mb_OrbBased_bodyTypeOfOrb(j2)

        if (i2bt .ne. j2bt) cycle

        orbTmps(:, :) = 0.0_R64
        call Method_Mb_OrbBased_ApplyInteractionOp(orbTmps, orbs, i2, j2, time)

        do j1 = 1, nOS
          do i1 = 1, nOS
            i1bt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)
            j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

            if (i1bt .ne. j1bt) cycle

            h2(i1, i2, j1, j2) = 0.5_R64 * Grid_InnerProduct(orbs(:, i1), orbTmps(:, j1))
          end do
        end do
      end do
    end do

    ! Handle the restricted case, if needed
    if (Orbs_restrictedQ) then
      h2(nOS + 1:, :nOS, nOS + 1:, :nOS) = h2(:nOS, :nOS, :nOS, :nOS)
      h2(:nOS, nOS + 1:, :nOS, nOS + 1:) = h2(:nOS, :nOS, :nOS, :nOS)
      h2(nOS + 1:, nOS + 1:, nOS + 1:, nOS + 1:) = h2(:nOS, :nOS, :nOS, :nOS)
    end if
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyKineticOp(dOrbs, orbs, time)
    use M_Method
    use M_Grid
    use M_Orbs
    use M_SysKinetic

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in)                  :: time

    complex(R64), allocatable :: orbTmp(:)

    integer(I32) :: nG, nOS
    integer(I32) :: i1, ibt

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    allocate (orbTmp(nG))

    do i1 = 1, nOS

      ibt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)

      orbTmp(:) = 0.0_R64
      call SysKinetic_MultiplyWithKineticOp(orbTmp, orbs(:, i1), time, ibt)

      dOrbs(:, i1) = dOrbs(:, i1) + orbTmp(:)

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyPotentialOp(dOrbs, orbs, time)
    use M_Grid
    use M_Grid_Ylm
    use M_Orbs
    use M_Method_Mb
    use M_SysPotential

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in) :: time

    integer(I32) :: nG, nOS
    integer(I32) :: iOrb, bt
    complex(R64), allocatable :: orbTmp(:)
    complex(R64), allocatable :: externalPotential(:)

    ! sizes
    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    allocate (orbTmp(nG))

    if (SysPotential_bodyTypeIndependentQ) then
      ! time-independent potential: compute once and apply to all orbitals
      call SysPotential_FillExternalPotential(externalPotential, time)
      do iOrb = 1, nOS
        call SysPotential_MultiplyWithExternalPotential(orbTmp, externalPotential, orbs(:, iOrb))
        dOrbs(:, iOrb) = dOrbs(:, iOrb) + orbTmp(:)
      end do
    else
      ! loop once per destination body type and apply to all orbitals of that type
      do bt = 1, Method_Mb_nBodyTypes
        call SysPotential_FillExternalPotential(externalPotential, time, bt_=bt)
        do iOrb = Method_Mb_OrbBased_nOrbsStart(bt), Method_Mb_OrbBased_nOrbsEnd(bt)
          call SysPotential_MultiplyWithExternalPotential(orbTmp, externalPotential, orbs(:, iOrb))
          dOrbs(:, iOrb) = dOrbs(:, iOrb) + orbTmp(:)
        end do
      end do
    end if

    deallocate (orbTmp)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyInteractionOp(dOrbs, orbs, i2, j2, time)
    use M_Grid
    use M_Grid_Ylm
    use M_Method_Mb
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    integer(I32), intent(in)               :: i2
    integer(I32), intent(in)               :: j2
    real(R64), intent(in)                  :: time

    integer(I32) :: nG, nOS
    integer(I32) :: iOrb, bt1, bt2
    complex(R64), allocatable :: src(:), interactionPotential(:), orbTmp(:)

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    allocate (orbTmp(nG))

    call SysInteraction_FillInteractionSrc(src, orbs(:, i2), orbs(:, j2))

    bt2 = Method_Mb_OrbBased_bodyTypeOfOrb(i2)

    if (SysInteraction_bodyTypeIndependentQ) then
      ! Compute once (bt1 and bt2 argument ignored by implementation) and apply to all orbitals
      call SysInteraction_FillInteractionPotential(interactionPotential, src, time, -1, -1)
      do iOrb = 1, nOS
        call SysInteraction_MultiplyWithInteractionPotential(orbTmp, interactionPotential, orbs(:, iOrb))
        dOrbs(:, iOrb) = dOrbs(:, iOrb) + orbTmp(:)
      end do
    else
      ! loop once per destination body type and apply to all orbitals of that type
      do bt1 = 1, Method_Mb_nBodyTypes
        call SysInteraction_FillInteractionPotential(interactionPotential, src, time, bt1, bt2)
        do iOrb = Method_Mb_OrbBased_nOrbsStart(bt1), Method_Mb_OrbBased_nOrbsEnd(bt1)
          call SysInteraction_MultiplyWithInteractionPotential(orbTmp, interactionPotential, orbs(:, iOrb))
          dOrbs(:, iOrb) = dOrbs(:, iOrb) + orbTmp(:)
        end do
      end do
    end if

    deallocate (src, interactionPotential, orbTmp)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyCorrelationOp(dOrbs, orbs, rdm1, rdm2, time, h2_)
    use M_Utils_RdmDiagonalize
    use M_Method
    use M_Grid
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    complex(R64), intent(in), contiguous             :: rdm1(:, :)
    complex(R64), intent(in), contiguous             :: rdm2(:, :, :, :)
    real(R64), intent(in)                  :: time
    complex(R64), intent(out), allocatable, optional :: h2_(:, :, :, :)

    real(R64), allocatable :: natocc(:)
    complex(R64), allocatable :: natorb(:, :)
    complex(R64), allocatable :: correlOp(:, :, :, :)

    complex(R64), allocatable :: orbTmps(:, :)

    complex(R64), allocatable :: rdm1spatial(:, :)
    complex(R64), allocatable :: rdm2spatial(:, :, :, :)

    integer(I32) :: nO, nOS
    integer(I32) :: j1, k1, k2, l2, i1
    integer(I32) :: k2bt, k1bt, j1bt, l2bt, i1bt
    real(R64) :: reg

    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    if (present(h2_)) then
      if (.not. allocated(h2_)) allocate (h2_(nO, nO, nO, nO))
      h2_(:, :, :, :) = 0.0_R64
    end if

    allocate (rdm1spatial(nOS, nOS))
    allocate (rdm2spatial(nOS, nOS, nOS, nOS))

    if (nO .eq. nOS) then
      rdm1spatial(:, :) = rdm1(1:nOS, 1:nOS)
      rdm2spatial(:, :, :, :) = rdm2(1:nOS, 1:nOS, 1:nOS, 1:nOS)

    else if (nO .eq. 2 * nOS) then
      rdm1spatial(:, :) = rdm1(1:nOS, 1:nOS) + rdm1(nOS + 1:, nOS + 1:)
      rdm2spatial(:, :, :, :) = rdm2(1:nOS, 1:nOS, 1:nOS, 1:nOS) + &
                                rdm2(nOS + 1:, :nOS, nOS + 1:, :nOS) + &
                                rdm2(:nOS, nOS + 1:, :nOS, nOS + 1:) + &
                                rdm2(nOS + 1:, nOS + 1:, nOS + 1:, nOS + 1:)
    else
      error stop "invalid nO to nOS relation in ApplyCorrelationOp"
    end if

    call RdmDiagonalize_Rdm1(natocc, natorb, rdm1spatial)

    reg = Method_Mb_OrbBased_regularizationParameter

    allocate (correlOp, mold=rdm2spatial)
    allocate (orbTmps, mold=orbs)

    correlOp(:, :, :, :) = 0.0_R64
    do j1 = 1, nOS
      do k1 = 1, nOS
        correlOp(:, :, j1, :) = correlOp(:, :, j1, :) + rdm2spatial(:, :, k1, :) * natorb(k1, j1)
      end do
    end do

    do j1 = 1, nOS
      correlOp(:, :, j1, :) = correlOp(:, :, j1, :) * natocc(j1) / (natocc(j1)**2 + reg)
    end do

    rdm2spatial(:, :, :, :) = 0.0_R64
    do k1 = 1, nOS
      do j1 = 1, nOS
        rdm2spatial(:, :, j1, :) = rdm2spatial(:, :, j1, :) + correlOp(:, :, k1, :) * conjg(natorb(j1, k1))
      end do
    end do

    do k2 = 1, nOS
      do l2 = 1, nOS
        l2bt = Method_Mb_OrbBased_bodyTypeOfOrb(l2)
        k2bt = Method_Mb_OrbBased_bodyTypeOfOrb(k2)

        if (l2bt .ne. k2bt) cycle

        orbTmps(:, :) = 0.0_R64
        call Method_Mb_OrbBased_ApplyInteractionOp(orbTmps, orbs, l2, k2, time)

        do j1 = 1, nOS
          j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

          do k1 = 1, nOS
            k1bt = Method_Mb_OrbBased_bodyTypeOfOrb(k1)

            if (j1bt .ne. k1bt) cycle

            dOrbs(:, j1) = dOrbs(:, j1) + orbTmps(:, k1) * rdm2spatial(k1, k2, j1, l2)
          end do
        end do

        if (present(h2_)) then
          do j1 = 1, nOS
            do i1 = 1, nOS
              i1bt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)
              j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

              if (i1bt .ne. j1bt) cycle

              h2_(i1, l2, j1, k2) = 0.5_R64 * Grid_InnerProduct(orbs(:, i1), orbTmps(:, j1))
            end do
          end do
        end if

      end do
    end do

    deallocate (rdm1spatial)
    deallocate (rdm2spatial)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyHartreeFockOp(dOrbs, orbs, rdm1, time)
    use M_Method
    use M_Grid
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    complex(R64), intent(in), contiguous             :: rdm1(:, :)
    real(R64), intent(in)                  :: time

    complex(R64), allocatable :: orbTmps(:, :)
    integer(I32) :: nOS
    integer(I32) :: j1, k2, l2
    integer(I32) :: j1bt, l2bt, k2bt

    allocate (orbTmps, mold=orbs)

    nOS = Orbs_nOrbsInState

    do k2 = 1, nOS
      do l2 = 1, nOS
        l2bt = Method_Mb_OrbBased_bodyTypeOfOrb(l2)
        k2bt = Method_Mb_OrbBased_bodyTypeOfOrb(k2)

        if (l2bt .ne. k2bt) cycle

        orbTmps(:, :) = 0.0_R64
        call Method_Mb_OrbBased_ApplyInteractionOp(orbTmps, orbs, l2, k2, time)

        do j1 = 1, nOS
          j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

          dOrbs(:, j1) = dOrbs(:, j1) + orbTmps(:, j1) * rdm1(k2, l2)

          if (k2bt .ne. j1bt) cycle
          dOrbs(:, k2) = dOrbs(:, k2) - orbTmps(:, j1) * rdm1(j1, l2)

        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplySingleBodyOp(dOrbs, orbs, time, h1_)
    use M_Method
    use M_Orbs
    use M_Grid
    use M_SysKinetic
    use M_SysPotential

    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    real(R64), intent(in)               :: time
    complex(R64), intent(out), allocatable, optional :: h1_(:, :)

    integer(I32) :: i1, j1, i1bt, j1bt, nOS, nO

    ! Get dimensions
    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    call Method_Mb_OrbBased_ApplyKineticOp(dOrbs, orbs, time)
    call Method_Mb_OrbBased_ApplyPotentialOp(dOrbs, orbs, time)

    if (present(h1_)) then

      if (.not. allocated(h1_)) allocate (h1_(nO, nO))
      h1_(:, :) = 0.0_R64

      do j1 = 1, nOS
        do i1 = 1, nOS

          i1bt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)
          j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

          if (i1bt .ne. j1bt) cycle

          h1_(i1, j1) = Grid_InnerProduct(orbs(:, i1), dOrbs(:, j1))

        end do
      end do

      if (Orbs_restrictedQ) then
        h1_(nOS + 1:, nOS + 1:) = h1_(:nOS, :nOS)
      end if

    end if

  end subroutine

end submodule
