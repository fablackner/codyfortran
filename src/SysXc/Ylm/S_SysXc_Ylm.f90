! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Ylm family of SysXc: angular product quadrature and its maps.
!>
!> @details
!> Tables (Setup): Gauss-Legendre nodes x_q, Θ_lm(x_q) for 0 ≤ m ≤ l ≤
!> lmaxTab (m up to mmaxTab + 1 for the derivatives), the analysis weights
!> 2π w_q Θ_lm(x_q), dΘ_lm/dθ(x_q) and 2π w_q dΘ_lm/dθ(x_q) (gradients),
!> the phases e^{imφ_k} and, on grid.ylm.fedvr, the band of the assembled
!> weak radial derivative B(i, j) = ∫ L_i L_j' dr.
!>
!> Transforms of one m (ChannelsToNodesM, NodesToChannelsM, with the table
!> to use): one real DGEMM on the complex arrays viewed as real,
!> (2 nRad × nL)(nL × nTheta) and back, with the sign (-1)^m of
!> Θ_l,-m = (-1)^m Θ_lm. Channel matrices of a
!> layout (orbitals, interaction source/potential) are gathered and scattered
!> with SphericalHarmonics_ChannelIndex. Loops over m (and over φ_k) are
!> OpenMP-parallel unless called in an active parallel region.
submodule(M_SysXc_Ylm) S_SysXc_Ylm
  use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer, c_associated
  use M_SysXc, only: SysXc_nPoints, SysXc_PARALLEL_MIN_WORK

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Radial points
  integer(I32) :: nRad = 0

  !> Largest l and |m| of the angular tables
  integer(I32) :: lmaxTab = 0, mmaxTab = 0

  !> Θ_lm(x_q), (0:lmaxTab, nTheta, 0:mmaxTab + 1), zero for l < m
  real(R64), allocatable, target :: thetaTab(:, :, :)
  !> 2π w_q Θ_lm(x_q), (nTheta, 0:lmaxTab, 0:mmaxTab)
  real(R64), allocatable, target :: wThetaTab(:, :, :)
  !> dΘ_lm/dθ (x_q), (0:lmaxTab, nTheta, 0:mmaxTab)
  real(R64), allocatable, target :: dThetaTab(:, :, :)
  !> 2π w_q dΘ_lm/dθ (x_q), (nTheta, 0:lmaxTab, 0:mmaxTab)
  real(R64), allocatable, target :: wdThetaTab(:, :, :)
  !> sin θ_q
  real(R64), allocatable :: sinTheta(:)
  !> e^{imφ_k}, (nPhi, -mmaxTab:mmaxTab)
  complex(R64), allocatable :: phases(:, :)

  !> Radial gradient (grid.ylm.fedvr): points r_i, FEDVR weights w_i (without
  !> r²), bands B(i, i + d) and B(i + d, i), d = -radialP..radialP
  real(R64), allocatable :: rPoints(:), radialDofWeights(:)
  real(R64), allocatable :: radialBand(:, :), radialAdjointBand(:, :)
  integer(I32) :: radialP = 0

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysXc_Ylm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm
    use M_SysInteraction_Ylm
    use M_SysXc
    use M_SysXc_Libxc

    integer(I32) :: nPhiDefault
    logical :: definiteMQ

    call Say_Fabricate("sysXc.ylm")

    if (.not. Json_GetExistence("grid.ylm")) error stop "grid.ylm is required for sysXc.ylm"
    if (Json_GetExistence("grid.ylm.fedvrEcs")) &
      error stop "sysXc.ylm: densities need a Hermitian radial grid (not grid.ylm.fedvrEcs)"

    !------------------------------------
    ! layout of the Hartree source and potential
    !------------------------------------

    SysXc_Ylm_lmaxSrc = 0
    SysXc_Ylm_mmaxSrc = 0
    if (SysXc_hartreeQ) then
      if (.not. Json_GetExistence("sysInteraction.ylm")) &
        error stop "sysXc.ylm: the Hartree potential requires sysInteraction.ylm (fabricate SysInteraction first)"
      SysXc_Ylm_lmaxSrc = SysInteraction_Ylm_lmax
      SysXc_Ylm_mmaxSrc = SysInteraction_Ylm_mmax
    end if

    !------------------------------------
    ! angular quadrature
    !------------------------------------

    definiteMQ = Json_GetExistence("method.mb.orbBased.mOfOrbs")
    nPhiDefault = 4 * Grid_Ylm_mmax + 1
    if (Grid_Ylm_mmax .eq. 0 .or. definiteMQ) nPhiDefault = 1

    SysXc_Ylm_nTheta = Json_Get("sysXc.ylm.nTheta", 2 * Grid_Ylm_lmax + 3)
    SysXc_Ylm_nPhi = Json_Get("sysXc.ylm.nPhi", nPhiDefault)
    if (SysXc_Ylm_nTheta < Grid_Ylm_lmax + 1) error stop "sysXc.ylm.nTheta must be at least grid.ylm.lmax + 1"
    if (SysXc_Ylm_nPhi < 1) error stop "sysXc.ylm.nPhi must be positive"

    SysXc_Ylm_axialQ = SysXc_Ylm_nPhi .eq. 1
    SysXc_Ylm_lmaxDensity = 2 * Grid_Ylm_lmax
    if (SysXc_Ylm_axialQ) then
      SysXc_Ylm_mmaxDensity = 0
      if (Grid_Ylm_mmax > 0 .and. .not. definiteMQ) &
        write (*, '(1X, A)') "note: sysXc.ylm.nPhi = 1 without method.mb.orbBased.mOfOrbs: axial approximation "// &
        "(phi-averaged density and an axially symmetric Kohn-Sham potential)"
    else
      SysXc_Ylm_mmaxDensity = 2 * Grid_Ylm_mmax
      if (SysXc_Ylm_nPhi < 2 * SysXc_Ylm_mmaxDensity + 1) &
        error stop "sysXc.ylm.nPhi must be 1 (axial mode) or at least 4 grid.ylm.mmax + 1"
      if (SysXc_Ylm_nPhi < 2 * SysXc_Ylm_mmaxSrc + 1) &
        error stop "sysXc.ylm.nPhi must be at least 2 sysInteraction.ylm.mmax + 1"
    end if

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysXc_dimension = 3
    SysXc_nPoints = Grid_Ylm_nRadial * SysXc_Ylm_nTheta * SysXc_Ylm_nPhi
    SysXc_Setup => Setup
    SysXc_AddDensity => AddDensity
    SysXc_MultiplyWithPotential => MultiplyWithPotential
    SysXc_FillInteractionSrc => FillInteractionSrc
    SysXc_AddInteractionPotential => AddInteractionPotential
    if (Json_GetExistence("grid.ylm.fedvr")) then
      SysXc_ApplyGradient => ApplyGradient
      SysXc_AddGradientAdjoint => AddGradientAdjoint
    end if

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysXc.ylm.libxc")) then
      call SysXc_Libxc_Fabricate("sysXc.ylm.libxc")

    else
      error stop "sysXc.ylm is missing one of: libxc"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Angular tables, weights of the xc points and (fedvr) the radial bands.
  subroutine Setup
    use stdlib_quadrature, only: gauss_legendre
    use M_Utils_Constants
    use M_Utils_Say
    use M_Grid_Ylm
    use M_SysXc

    real(R64), allocatable :: x(:), w(:)
    integer(I32) :: nT, nPh, l, m, q, k, iRad, p

    call Say_Setup("sysXc.ylm")

    if (.not. allocated(Grid_Ylm_radialWeights)) error stop "sysXc.ylm: call Grid_Setup before SysXc_Setup"

    nRad = Grid_Ylm_nRadial
    nT = SysXc_Ylm_nTheta
    nPh = SysXc_Ylm_nPhi
    lmaxTab = max(Grid_Ylm_lmax, SysXc_Ylm_lmaxSrc, SysXc_Ylm_lmaxDensity)
    mmaxTab = max(Grid_Ylm_mmax, SysXc_Ylm_mmaxSrc, SysXc_Ylm_mmaxDensity)
    mmaxTab = min(mmaxTab, lmaxTab)

    !------------------------------------
    ! Gauss-Legendre nodes, exactly symmetric
    !------------------------------------

    allocate (x(nT), w(nT))
    call gauss_legendre(x, w)
    x(:) = 0.5_R64 * (x(:) - x(nT:1:-1))
    w(:) = 0.5_R64 * (w(:) + w(nT:1:-1))
    sinTheta = sqrt(1.0_R64 - x(:)**2)
    SysXc_Ylm_cosTheta = x

    !------------------------------------
    ! normalized associated Legendre functions (recurrences of Grid_Ylm):
    ! Θ_00 = 1/sqrt(4π), Θ_mm = -sqrt((2m+1)/(2m)) sinθ Θ_m-1,m-1,
    ! Θ_m+1,m = sqrt(2m+3) x Θ_mm, Θ_lm = a_lm (x Θ_l-1,m - b_lm Θ_l-2,m)
    !------------------------------------

    if (allocated(thetaTab)) deallocate (thetaTab, wThetaTab, dThetaTab, wdThetaTab)
    allocate (thetaTab(0:lmaxTab, nT, 0:mmaxTab + 1))
    thetaTab(:, :, :) = 0.0_R64
    thetaTab(0, :, 0) = 1.0_R64 / sqrt(4.0_R64 * PI)
    do m = 1, min(mmaxTab + 1, lmaxTab)
      thetaTab(m, :, m) = -sqrt((2 * m + 1) / (2.0_R64 * m)) * sinTheta(:) * thetaTab(m - 1, :, m - 1)
    end do
    do m = 0, min(mmaxTab + 1, lmaxTab)
      if (m + 1 <= lmaxTab) thetaTab(m + 1, :, m) = sqrt(2.0_R64 * m + 3.0_R64) * x(:) * thetaTab(m, :, m)
      do l = m + 2, lmaxTab
        thetaTab(l, :, m) = sqrt((4.0_R64 * l**2 - 1.0_R64) / real(l**2 - m**2, kind=R64)) * &
                            (x(:) * thetaTab(l - 1, :, m) - &
                             sqrt(real((l - 1)**2 - m**2, kind=R64) / (4.0_R64 * (l - 1)**2 - 1.0_R64)) * &
                             thetaTab(l - 2, :, m))
      end do
    end do

    ! θ derivatives from the ladder relations (Θ_l,-1 = -Θ_l1):
    ! dΘ_lm/dθ = ½ [sqrt((l-m)(l+m+1)) Θ_l,m+1 - sqrt((l+m)(l-m+1)) Θ_l,m-1],
    ! dΘ_l0/dθ = sqrt(l(l+1)) Θ_l1
    allocate (dThetaTab(0:lmaxTab, nT, 0:mmaxTab))
    dThetaTab(:, :, :) = 0.0_R64
    do m = 0, mmaxTab
      do l = m, lmaxTab
        if (m .eq. 0) then
          if (l > 0) dThetaTab(l, :, 0) = sqrt(real(l * (l + 1), kind=R64)) * thetaTab(l, :, 1)
        else
          dThetaTab(l, :, m) = 0.5_R64 * (sqrt(real((l - m) * (l + m + 1), kind=R64)) * thetaTab(l, :, m + 1) - &
                                          sqrt(real((l + m) * (l - m + 1), kind=R64)) * thetaTab(l, :, m - 1))
        end if
      end do
    end do

    allocate (wThetaTab(nT, 0:lmaxTab, 0:mmaxTab), wdThetaTab(nT, 0:lmaxTab, 0:mmaxTab))
    do m = 0, mmaxTab
      do l = 0, lmaxTab
        wThetaTab(:, l, m) = 2.0_R64 * PI * w(:) * thetaTab(l, :, m)
        wdThetaTab(:, l, m) = 2.0_R64 * PI * w(:) * dThetaTab(l, :, m)
      end do
    end do

    !------------------------------------
    ! phases and weights of the xc points
    !------------------------------------

    if (allocated(phases)) deallocate (phases)
    allocate (phases(nPh, -mmaxTab:mmaxTab))
    do m = -mmaxTab, mmaxTab
      do k = 1, nPh
        phases(k, m) = exp(cmplx(0.0_R64, 2.0_R64 * PI * m * (k - 1) / real(nPh, kind=R64), kind=R64))
      end do
    end do

    SysXc_nPoints = nRad * nT * nPh
    if (allocated(SysXc_weights)) deallocate (SysXc_weights)
    allocate (SysXc_weights(SysXc_nPoints))
    do k = 1, nPh
      do q = 1, nT
        do iRad = 1, nRad
          p = iRad + nRad * (q - 1) + nRad * nT * (k - 1)
          SysXc_weights(p) = Grid_Ylm_radialWeights(iRad) * 2.0_R64 * PI * w(q) / real(nPh, kind=R64)
        end do
      end do
    end do

    if (associated(SysXc_ApplyGradient)) call SetupRadialBands

    write (*, '(1X, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, 2A)') "sysXc.ylm: ", SysXc_nPoints, " xc points (", nRad, &
      " radial x ", nT, " theta x ", nPh, " phi), tables l <= ", lmaxTab, ", |m| <= ", mmaxTab, &
      ", mode ", merge("axial  ", "general", SysXc_Ylm_axialQ)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Band of the assembled weak radial derivative B(i, j) = Σ_e wRef(i)
  !> dL_j/dξ(ξ_i) on the radial points of grid.ylm.fedvr (r = 0 excluded: its
  !> row and column drop out, r f vanishes there).
  subroutine SetupRadialBands
    use M_Grid_Ylm
    use M_Grid_Ylm_Fedvr, only: fedvrCtx => Grid_Ylm_Fedvr_fedvrCtx, derivativeCtx => Grid_Ylm_Fedvr_derivativeCtx

    real(R64), allocatable :: weak(:, :)
    integer(I32) :: nLoc, nE, iE, iLocalStart, iLocalEnd, offset, iLoc, jLoc, i, d

    nLoc = fedvrCtx % nLocals
    nE = fedvrCtx % nElements
    radialP = nLoc - 1
    if (fedvrCtx % nPoints .ne. nRad) error stop "sysXc.ylm: radial FEDVR points differ from Grid_Ylm_nRadial"

    allocate (weak(nRad, nRad))
    weak(:, :) = 0.0_R64
    do iE = 1, nE
      associate (element => fedvrCtx % elements(iE))
        iLocalStart = 1
        iLocalEnd = nLoc
        if (fedvrCtx % xminExcludedQ .and. iE .eq. 1) iLocalStart = 2
        if (fedvrCtx % xmaxExcludedQ .and. iE .eq. nE) iLocalEnd = nLoc - 1
        offset = element % iStart - iLocalStart
        do jLoc = iLocalStart, iLocalEnd
          do iLoc = iLocalStart, iLocalEnd
            weak(offset + iLoc, offset + jLoc) = weak(offset + iLoc, offset + jLoc) + &
                                                 derivativeCtx % firstOrderMatrixWeak(iLoc, jLoc)
          end do
        end do
      end associate
    end do

    rPoints = Grid_Ylm_radialPoints
    radialDofWeights = fedvrCtx % weights
    if (allocated(radialBand)) deallocate (radialBand, radialAdjointBand)
    allocate (radialBand(-radialP:radialP, nRad), radialAdjointBand(-radialP:radialP, nRad))
    radialBand(:, :) = 0.0_R64
    radialAdjointBand(:, :) = 0.0_R64
    do i = 1, nRad
      do d = max(-radialP, 1 - i), min(radialP, nRad - i)
        radialBand(d, i) = weak(i, i + d)
        radialAdjointBand(d, i) = weak(i + d, i)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> density += Re Σ_l conj(φ_l) g_l on the xc points.
  subroutine AddDensity(density, orbs, weightedOrbs)
    use M_Grid_Ylm

    real(R64), intent(inout), contiguous, target :: density(:)
    complex(R64), intent(in), contiguous, target :: orbs(:, :)
    complex(R64), intent(in), contiguous, target :: weightedOrbs(:, :)

    complex(R64), allocatable, target :: valsA(:, :, :), valsB(:, :, :), fM(:, :), gA(:, :), gB(:, :)
    real(R64), allocatable :: densityLoc(:)
    integer(I32) :: nT, nPh, nO, l, m, mm, iPair, nPairs, nL, nNode
    logical :: sameQ

    call CheckReady("SysXc_AddDensity")
    nT = SysXc_Ylm_nTheta
    nPh = SysXc_Ylm_nPhi
    nO = size(orbs, 2)
    nNode = nRad * nT
    sameQ = c_associated(c_loc(orbs), c_loc(weightedOrbs))
    if (size(density) .ne. nNode * nPh) error stop "SysXc_AddDensity (ylm): density must have SysXc_nPoints entries"

    if (SysXc_Ylm_axialQ) then
      ! Every (orbital, m) block separately: density += Re conj(f_m) g_m
      nPairs = nO * (2 * Grid_Ylm_mmax + 1)
      allocate (fM(nRad, Grid_Ylm_lmax + 1), gA(nRad, nT), gB(nRad, nT), densityLoc(nNode))
      densityLoc(:) = 0.0_R64
      do iPair = 1, nPairs
        l = (iPair - 1) / (2 * Grid_Ylm_mmax + 1) + 1
        m = mod(iPair - 1, 2 * Grid_Ylm_mmax + 1) - Grid_Ylm_mmax
        mm = abs(m)
        nL = Grid_Ylm_lmax - mm + 1
        call GatherM(fM(:, 1:nL), orbs(:, l), m, Grid_Ylm_mmax)
        call ChannelsToNodesM(gA, fM(:, 1:nL), m, thetaTab)
        if (.not. sameQ) then
          call GatherM(fM(:, 1:nL), weightedOrbs(:, l), m, Grid_Ylm_mmax)
          call ChannelsToNodesM(gB, fM(:, 1:nL), m, thetaTab)
        end if
        call AccumulateProduct(densityLoc, gA, gB, sameQ)
      end do
      density(1:nNode) = density(1:nNode) + densityLoc(:)
      return
    end if

    ! General mode: orbitals on all points
    allocate (valsA(nRad, nT, nPh), valsB(nRad, nT, nPh))
    do l = 1, nO
      call OrbitalToPoints(valsA, orbs(:, l))
      if (.not. sameQ) call OrbitalToPoints(valsB, weightedOrbs(:, l))
      call AccumulateProduct(density, valsA, valsB, sameQ)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> density += Re conj(a) b (or |a|² if sameQ) point by point (a, b of
  !> rank 2 or 3 in the layout of the xc points; b unused if sameQ).
  subroutine AccumulateProduct(density, a, b, sameQ)
!$  use omp_lib, only: omp_in_parallel

    real(R64), intent(inout), contiguous :: density(:)
    complex(R64), intent(in), contiguous, target :: a(..)
    complex(R64), intent(in), contiguous, target :: b(..)
    logical, intent(in) :: sameQ

    complex(R64), pointer :: a1(:), b1(:)
    integer(I32) :: p

    if (size(a) .ne. size(density) .or. size(b) .ne. size(density)) &
      error stop "sysXc.ylm: AccumulateProduct sizes differ"
    call c_f_pointer(c_loc(a), a1, [size(a)])
    call c_f_pointer(c_loc(b), b1, [size(b)])

    !$omp parallel do default(shared) private(p) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do p = 1, size(density)
      if (sameQ) then
        density(p) = density(p) + real(a1(p), kind=R64)**2 + aimag(a1(p))**2
      else
        density(p) = density(p) + real(a1(p), kind=R64) * real(b1(p), kind=R64) + aimag(a1(p)) * aimag(b1(p))
      end if
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> resOrb = v orb: nodes, multiplication, projection (per m block in the
  !> axial mode).
  subroutine MultiplyWithPotential(resOrb, pot, orb)
!$  use omp_lib, only: omp_in_parallel
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous, target :: resOrb(:)
    real(R64), intent(in), contiguous, target :: pot(:)
    complex(R64), intent(in), contiguous, target :: orb(:)

    complex(R64), allocatable, target :: vals(:, :, :), fM(:, :), g(:, :)
    real(R64), contiguous, pointer :: pot2(:, :), pot3(:, :, :)
    integer(I32) :: nT, nPh, m, nL, k

    call CheckReady("SysXc_MultiplyWithPotential")
    nT = SysXc_Ylm_nTheta
    nPh = SysXc_Ylm_nPhi

    if (SysXc_Ylm_axialQ) then
      pot2(1:nRad, 1:nT) => pot
      do m = -Grid_Ylm_mmax, Grid_Ylm_mmax
        nL = Grid_Ylm_lmax - abs(m) + 1
        allocate (fM(nRad, nL), g(nRad, nT))
        call GatherM(fM, orb, m, Grid_Ylm_mmax)
        call ChannelsToNodesM(g, fM, m, thetaTab)
        g(:, :) = pot2(:, :) * g(:, :)
        call NodesToChannelsM(fM, g, m, wThetaTab)
        call ScatterM(resOrb, fM, m, Grid_Ylm_mmax)
        deallocate (fM, g)
      end do
      return
    end if

    allocate (vals(nRad, nT, nPh))
    call OrbitalToPoints(vals, orb)
    pot3(1:nRad, 1:nT, 1:nPh) => pot
    !$omp parallel do default(shared) private(k) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do k = 1, nPh
      vals(:, :, k) = pot3(:, :, k) * vals(:, :, k)
    end do
    !$omp end parallel do
    call PointsToOrbital(resOrb, vals)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> vals(:, :, k) = Σ_m f_m(r, θ_q) e^{imφ_k} of an orbital (general mode).
  subroutine OrbitalToPoints(vals, orb)
!$  use omp_lib, only: omp_in_parallel
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous :: vals(:, :, :)
    complex(R64), intent(in), contiguous :: orb(:)

    complex(R64), allocatable :: gAll(:, :, :), fM(:, :)
    integer(I32) :: mmax, m, nL, k

    mmax = Grid_Ylm_mmax
    allocate (gAll(nRad, SysXc_Ylm_nTheta, -mmax:mmax))

    do m = -mmax, mmax
      nL = Grid_Ylm_lmax - abs(m) + 1
      allocate (fM(nRad, nL))
      call GatherM(fM, orb, m, mmax)
      call ChannelsToNodesM(gAll(:, :, m), fM, m, thetaTab)
      deallocate (fM)
    end do

    !$omp parallel do default(shared) private(k, m) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do k = 1, SysXc_Ylm_nPhi
      vals(:, :, k) = gAll(:, :, 0)
      do m = 1, mmax
        vals(:, :, k) = vals(:, :, k) + gAll(:, :, m) * phases(k, m) + gAll(:, :, -m) * phases(k, -m)
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> orb_lm = Σ_q,k (2π w_q/nPhi) Θ_lm(x_q) e^{-imφ_k} vals (general mode).
  subroutine PointsToOrbital(orb, vals)
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous :: orb(:)
    complex(R64), intent(in), contiguous :: vals(:, :, :)

    complex(R64), allocatable :: uM(:, :), fM(:, :)
    integer(I32) :: m, nL, k

    do m = -Grid_Ylm_mmax, Grid_Ylm_mmax
      nL = Grid_Ylm_lmax - abs(m) + 1
      allocate (uM(nRad, SysXc_Ylm_nTheta), fM(nRad, nL))
      uM(:, :) = vals(:, :, 1) * conjg(phases(1, m))
      do k = 2, SysXc_Ylm_nPhi
        uM(:, :) = uM(:, :) + vals(:, :, k) * conjg(phases(k, m))
      end do
      uM(:, :) = uM(:, :) / real(SysXc_Ylm_nPhi, kind=R64)
      call NodesToChannelsM(fM, uM, m, wThetaTab)
      call ScatterM(orb, fM, m, Grid_Ylm_mmax)
      deallocate (uM, fM)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Interaction source of a density: channels l ≤ lmaxSrc, |m| ≤ mmaxSrc
  !> (m = 0 in the axial mode) of the projection onto Y_lm, times the radial
  !> metric weights (the convention of SysInteraction_FillInteractionSrc).
  subroutine FillInteractionSrc(src, density)
    use M_Grid_Ylm
    use M_SysInteraction_Ylm

    complex(R64), intent(out), allocatable, target :: src(:)
    real(R64), intent(in), contiguous, target :: density(:)

    real(R64), contiguous, pointer :: density3(:, :, :)
    complex(R64), allocatable :: rhoAll(:, :, :), cM(:, :)
    integer(I32) :: m, mMax, nL, iL

    call CheckReady("SysXc_FillInteractionSrc")
    allocate (src(SysInteraction_Ylm_nPoints))
    src(:) = (0.0_R64, 0.0_R64)
    density3(1:nRad, 1:SysXc_Ylm_nTheta, 1:SysXc_Ylm_nPhi) => density

    mMax = SysXc_Ylm_mmaxSrc
    if (SysXc_Ylm_axialQ) mMax = 0

    ! Fourier coefficients of the real density: rho_-m = conj(rho_m)
    allocate (rhoAll(nRad, SysXc_Ylm_nTheta, 0:mMax))
    call FillFourierCoefficients(rhoAll, density3, mMax)
    do m = -mMax, mMax
      nL = SysXc_Ylm_lmaxSrc - abs(m) + 1
      allocate (cM(nRad, nL))
      if (m < 0) then
        call NodesToChannelsM(cM, conjg(rhoAll(:, :, -m)), m, wThetaTab)
      else
        call NodesToChannelsM(cM, rhoAll(:, :, m), m, wThetaTab)
      end if
      do iL = 1, nL
        cM(:, iL) = cM(:, iL) * Grid_Ylm_radialMetricWeights(:)
      end do
      call ScatterM(src, cM, m, SysXc_Ylm_mmaxSrc)
      deallocate (cM)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> pot += Re V on the xc points for V in the interaction layout (m = 0 in
  !> the axial mode).
  subroutine AddInteractionPotential(pot, interactionPotential)
!$  use omp_lib, only: omp_in_parallel

    real(R64), intent(inout), contiguous, target :: pot(:)
    complex(R64), intent(in), contiguous, target :: interactionPotential(:)

    real(R64), contiguous, pointer :: pot3(:, :, :)
    complex(R64), allocatable :: gAll(:, :, :), vM(:, :)
    integer(I32) :: m, mMax, nL, k

    call CheckReady("SysXc_AddInteractionPotential")
    pot3(1:nRad, 1:SysXc_Ylm_nTheta, 1:SysXc_Ylm_nPhi) => pot

    mMax = SysXc_Ylm_mmaxSrc
    if (SysXc_Ylm_axialQ) mMax = 0
    allocate (gAll(nRad, SysXc_Ylm_nTheta, -mMax:mMax))

    do m = -mMax, mMax
      nL = SysXc_Ylm_lmaxSrc - abs(m) + 1
      allocate (vM(nRad, nL))
      call GatherM(vM, interactionPotential, m, SysXc_Ylm_mmaxSrc)
      call ChannelsToNodesM(gAll(:, :, m), vM, m, thetaTab)
      deallocate (vM)
    end do

    !$omp parallel do default(shared) private(k, m) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do k = 1, SysXc_Ylm_nPhi
      pot3(:, :, k) = pot3(:, :, k) + real(gAll(:, :, 0), kind=R64)
      do m = 1, mMax
        pot3(:, :, k) = pot3(:, :, k) + real(gAll(:, :, m) * phases(k, m) + gAll(:, :, -m) * phases(k, -m), kind=R64)
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> gradF = (∂_r f, (1/r) ∂_θ f, (1/(r sinθ)) ∂_φ f) of a real field with
  !> angular band limit lmaxDensity, mmaxDensity (module header).
  subroutine ApplyGradient(gradF, f)
!$  use omp_lib, only: omp_in_parallel

    real(R64), intent(out), contiguous, target :: gradF(:, :)
    real(R64), intent(in), contiguous, target :: f(:)

    real(R64), contiguous, pointer :: f3(:, :, :), gTheta(:, :, :), gPhi(:, :, :)
    complex(R64), allocatable :: fAll(:, :, :), dAll(:, :, :), cM(:, :)
    integer(I32) :: nT, nPh, mD, m, nL, k, q, iCol

    call CheckReady("SysXc_ApplyGradient")
    nT = SysXc_Ylm_nTheta
    nPh = SysXc_Ylm_nPhi
    mD = SysXc_Ylm_mmaxDensity

    f3(1:nRad, 1:nT, 1:nPh) => f
    gTheta(1:nRad, 1:nT, 1:nPh) => gradF(:, 2)
    gPhi(1:nRad, 1:nT, 1:nPh) => gradF(:, 3)

    ! Radial component, column by column
    !$omp parallel do default(shared) private(iCol) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do iCol = 1, nT * nPh
      call RadialDerivative(gradF((iCol - 1) * nRad + 1:iCol * nRad, 1), f((iCol - 1) * nRad + 1:iCol * nRad))
    end do
    !$omp end parallel do

    ! Angular components from the Fourier coefficients f_m (m >= 0: f real)
    allocate (fAll(nRad, nT, 0:mD), dAll(nRad, nT, 0:mD))
    call FillFourierCoefficients(fAll, f3, mD)
    do m = 0, mD
      nL = SysXc_Ylm_lmaxDensity - m + 1
      allocate (cM(nRad, nL))
      call NodesToChannelsM(cM, fAll(:, :, m), m, wThetaTab)
      call ChannelsToNodesM(dAll(:, :, m), cM, m, dThetaTab)
      deallocate (cM)
    end do

    !$omp parallel do default(shared) private(k, m, q) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do k = 1, nPh
      gTheta(:, :, k) = real(dAll(:, :, 0), kind=R64)
      gPhi(:, :, k) = 0.0_R64
      do m = 1, mD
        gTheta(:, :, k) = gTheta(:, :, k) + 2.0_R64 * real(dAll(:, :, m) * phases(k, m), kind=R64)
        gPhi(:, :, k) = gPhi(:, :, k) - 2.0_R64 * m * aimag(fAll(:, :, m) * phases(k, m))
      end do
      do q = 1, nT
        gTheta(:, q, k) = gTheta(:, q, k) / rPoints(:)
        gPhi(:, q, k) = gPhi(:, q, k) / (rPoints(:) * sinTheta(q))
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> f += G^† gradF in the metric of the xc points: the transposed radial
  !> band, analysis with dΘ/dθ and synthesis with Θ (θ), and -∂_φ (φ).
  subroutine AddGradientAdjoint(f, gradF)
!$  use omp_lib, only: omp_in_parallel

    real(R64), intent(inout), contiguous, target :: f(:)
    real(R64), intent(in), contiguous, target :: gradF(:, :)

    real(R64), contiguous, pointer :: f3(:, :, :)
    real(R64), allocatable :: uTheta(:, :, :), uPhi(:, :, :)
    complex(R64), allocatable :: thetaAll(:, :, :), phiAll(:, :, :), backAll(:, :, :), cM(:, :)
    real(R64), allocatable :: column(:)
    integer(I32) :: nT, nPh, mD, m, nL, k, q, iCol

    call CheckReady("SysXc_AddGradientAdjoint")
    nT = SysXc_Ylm_nTheta
    nPh = SysXc_Ylm_nPhi
    mD = SysXc_Ylm_mmaxDensity
    f3(1:nRad, 1:nT, 1:nPh) => f

    ! Radial component
    !$omp parallel default(shared) private(iCol, column) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    allocate (column(nRad))
    !$omp do schedule(static)
    do iCol = 1, nT * nPh
      call RadialDerivativeAdjoint(column, gradF((iCol - 1) * nRad + 1:iCol * nRad, 1))
      f((iCol - 1) * nRad + 1:iCol * nRad) = f((iCol - 1) * nRad + 1:iCol * nRad) + column(:)
    end do
    !$omp end do
    deallocate (column)
    !$omp end parallel

    ! Angular components: u_θ/r and u_φ/(r sinθ), then their Fourier coefficients
    allocate (uTheta(nRad, nT, nPh), uPhi(nRad, nT, nPh))
    !$omp parallel do default(shared) private(k, q) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do k = 1, nPh
      do q = 1, nT
        uTheta(:, q, k) = gradF(nRad * (q - 1) + nRad * nT * (k - 1) + 1:nRad * q + nRad * nT * (k - 1), 2) / rPoints(:)
        uPhi(:, q, k) = gradF(nRad * (q - 1) + nRad * nT * (k - 1) + 1:nRad * q + nRad * nT * (k - 1), 3) / &
                        (rPoints(:) * sinTheta(q))
      end do
    end do
    !$omp end parallel do

    allocate (thetaAll(nRad, nT, 0:mD), phiAll(nRad, nT, 0:mD), backAll(nRad, nT, 0:mD))
    call FillFourierCoefficients(thetaAll, uTheta, mD)
    call FillFourierCoefficients(phiAll, uPhi, mD)

    do m = 0, mD
      nL = SysXc_Ylm_lmaxDensity - m + 1
      allocate (cM(nRad, nL))
      call NodesToChannelsM(cM, thetaAll(:, :, m), m, wdThetaTab)
      call ChannelsToNodesM(backAll(:, :, m), cM, m, thetaTab)
      deallocate (cM)
    end do

    !$omp parallel do default(shared) private(k, m) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do k = 1, nPh
      f3(:, :, k) = f3(:, :, k) + real(backAll(:, :, 0), kind=R64)
      do m = 1, mD
        f3(:, :, k) = f3(:, :, k) + 2.0_R64 * real(backAll(:, :, m) * phases(k, m), kind=R64) &
                      + 2.0_R64 * m * aimag(phiAll(:, :, m) * phases(k, m))
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> fAll(:, :, m) = (1/nPhi) Σ_k f(:, :, k) e^{-imφ_k}, m = 0..mD.
  subroutine FillFourierCoefficients(fAll, f3, mD)
!$  use omp_lib, only: omp_in_parallel

    complex(R64), intent(out), contiguous :: fAll(:, :, 0:)
    real(R64), intent(in), contiguous :: f3(:, :, :)
    integer(I32), intent(in) :: mD

    integer(I32) :: m, k

    !$omp parallel do default(shared) private(m, k) schedule(static) &
    !$omp   if (.not. omp_in_parallel() .and. SysXc_nPoints >= SysXc_PARALLEL_MIN_WORK)
    do m = 0, mD
      fAll(:, :, m) = f3(:, :, 1) * conjg(phases(1, m))
      do k = 2, size(f3, 3)
        fAll(:, :, m) = fAll(:, :, m) + f3(:, :, k) * conjg(phases(k, m))
      end do
      fAll(:, :, m) = fAll(:, :, m) / real(size(f3, 3), kind=R64)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> df = ∂_r f = (D(r f) - f)/r on one radial column, D = w^-1 B.
  pure subroutine RadialDerivative(df, f)
    real(R64), intent(out) :: df(:)
    real(R64), intent(in) :: f(:)

    integer(I32) :: i, d
    real(R64) :: acc

    do i = 1, nRad
      acc = 0.0_R64
      do d = max(-radialP, 1 - i), min(radialP, nRad - i)
        acc = acc + radialBand(d, i) * rPoints(i + d) * f(i + d)
      end do
      df(i) = (acc / radialDofWeights(i) - f(i)) / rPoints(i)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Adjoint of RadialDerivative in the radial metric w r²:
  !> (G^† u)_j = (1/(w_j r_j)) Σ_i B(i, j) r_i u_i - u_j / r_j.
  pure subroutine RadialDerivativeAdjoint(res, u)
    real(R64), intent(out) :: res(:)
    real(R64), intent(in) :: u(:)

    integer(I32) :: j, d
    real(R64) :: acc

    do j = 1, nRad
      acc = 0.0_R64
      do d = max(-radialP, 1 - j), min(radialP, nRad - j)
        acc = acc + radialAdjointBand(d, j) * rPoints(j + d) * u(j + d)
      end do
      res(j) = acc / (radialDofWeights(j) * rPoints(j)) - u(j) / rPoints(j)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> nodes(r, q) = Σ_iL fM(r, iL) tab(|m| + iL - 1, q, |m|), tab = thetaTab
  !> (Θ) or dThetaTab (dΘ/dθ), times (-1)^m for m < 0.
  subroutine ChannelsToNodesM(nodes, fM, m, tab)
    use M_Utils_BlasLib, only: DGEMM

    complex(R64), intent(out), contiguous, target :: nodes(:, :)
    complex(R64), intent(in), contiguous, target :: fM(:, :)
    integer(I32), intent(in) :: m
    real(R64), intent(in) :: tab(0:lmaxTab, SysXc_Ylm_nTheta, 0:*)

    real(R64), contiguous, pointer :: nodesR(:, :), fMR(:, :)
    real(R64) :: sgn
    integer(I32) :: nRad2, nL, mm, nT

    mm = abs(m)
    nL = size(fM, 2)
    nT = SysXc_Ylm_nTheta
    nRad2 = 2 * size(nodes, 1)
    if (nL .eq. 0) then
      nodes(:, :) = (0.0_R64, 0.0_R64)
      return
    end if
    if (mm + nL - 1 > lmaxTab .or. mm > mmaxTab) error stop "sysXc.ylm: channel outside the angular tables"

    sgn = merge(-1.0_R64, 1.0_R64, m < 0 .and. mod(mm, 2) .eq. 1)
    call c_f_pointer(c_loc(nodes), nodesR, [nRad2, nT])
    call c_f_pointer(c_loc(fM), fMR, [nRad2, nL])

    call DGEMM('N', 'N', nRad2, nT, nL, sgn, fMR, nRad2, tab(mm, 1, mm), lmaxTab + 1, 0.0_R64, nodesR, nRad2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> fM(r, iL) = Σ_q nodes(r, q) wTab(q, |m| + iL - 1, |m|), wTab = wThetaTab
  !> (2π w_q Θ) or wdThetaTab (2π w_q dΘ/dθ), times (-1)^m for m < 0.
  subroutine NodesToChannelsM(fM, nodes, m, wTab)
    use M_Utils_BlasLib, only: DGEMM

    complex(R64), intent(out), contiguous, target :: fM(:, :)
    complex(R64), intent(in), contiguous, target :: nodes(:, :)
    integer(I32), intent(in) :: m
    real(R64), intent(in) :: wTab(SysXc_Ylm_nTheta, 0:lmaxTab, 0:*)

    real(R64), contiguous, pointer :: nodesR(:, :), fMR(:, :)
    real(R64) :: sgn
    integer(I32) :: nRad2, nL, mm, nT

    mm = abs(m)
    nL = size(fM, 2)
    nT = SysXc_Ylm_nTheta
    nRad2 = 2 * size(nodes, 1)
    if (nL .eq. 0) return
    if (mm + nL - 1 > lmaxTab .or. mm > mmaxTab) error stop "sysXc.ylm: channel outside the angular tables"

    sgn = merge(-1.0_R64, 1.0_R64, m < 0 .and. mod(mm, 2) .eq. 1)
    call c_f_pointer(c_loc(nodes), nodesR, [nRad2, nT])
    call c_f_pointer(c_loc(fM), fMR, [nRad2, nL])

    call DGEMM('N', 'N', nRad2, nL, nT, sgn, nodesR, nRad2, wTab(1, mm, mm), nT, 0.0_R64, fMR, nRad2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> fM(:, iL) = channel (|m| + iL - 1, m) of f in a layout with |m| ≤ mmaxL.
  subroutine GatherM(fM, f, m, mmaxL)
    use M_Utils_SphericalHarmonics, only: SphericalHarmonics_ChannelIndex

    complex(R64), intent(out), contiguous :: fM(:, :)
    complex(R64), intent(in), contiguous :: f(:)
    integer(I32), intent(in) :: m
    integer(I32), intent(in) :: mmaxL

    integer(I32) :: iL, l, i0

    do iL = 1, size(fM, 2)
      l = abs(m) + iL - 1
      i0 = (SphericalHarmonics_ChannelIndex(l, m, mmaxL) - 1) * nRad
      fM(:, iL) = f(i0 + 1:i0 + nRad)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Channel (|m| + iL - 1, m) of f = fM(:, iL) in a layout with |m| ≤ mmaxL.
  subroutine ScatterM(f, fM, m, mmaxL)
    use M_Utils_SphericalHarmonics, only: SphericalHarmonics_ChannelIndex

    complex(R64), intent(inout), contiguous :: f(:)
    complex(R64), intent(in), contiguous :: fM(:, :)
    integer(I32), intent(in) :: m
    integer(I32), intent(in) :: mmaxL

    integer(I32) :: iL, l, i0

    do iL = 1, size(fM, 2)
      l = abs(m) + iL - 1
      i0 = (SphericalHarmonics_ChannelIndex(l, m, mmaxL) - 1) * nRad
      f(i0 + 1:i0 + nRad) = fM(:, iL)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine CheckReady(caller)
    character(len=*), intent(in) :: caller

    if (.not. allocated(thetaTab)) error stop caller//" (ylm): SysXc_Setup was not called"

  end subroutine

end submodule
