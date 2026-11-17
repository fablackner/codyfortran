! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Implementation submodule of the Hubbard coefficient backend.
!>
!> - `Coeffs_Hubbard_Fabricate`: checks, configuration counts, pointer binding
!>   and dispatch to the spin-symmetry variant (which sets the parity and the
!>   number of coefficients)
!> - `Setup`: configurations of each spin and their single excitations
!> - `ApplyH1FillRdm1`, `ApplyH2FillRdm2`: Hamiltonian application and RDMs on
!>   the product space c(iUp, iDn); a sector vector is expanded to the product
!>   space and the result projected back
!> - `ApplyExcitation`: ladder operators of one spin (product space only)
!> - the index maps of the three variants
!>
!> RDM conventions are those of coeffs.generic (M_Coeffs): rdm1(p,q) =
!> <a+_q a_p>, rdm2(p,q,r,s) = <a+_r a+_s a_q a_p> over the spin orbitals
!> (up orbitals first), the two spins treated as one fermion species
!> (antisymmetric also under the exchange of an up and a down index).
submodule(M_Coeffs_Hubbard) S_Coeffs_Hubbard
  implicit none

  !> binomials(n, k) = binom(n, k) for 0 <= n, k <= nSites (zero for k > n),
  !> the colex ranks of the configurations.
  integer(I64), allocatable :: binomials(:, :)

  !> H1 of one spin for the h1 of one call: diag(i) = <i|H1|i> and the
  !> nonzero off-diagonal elements <i|H1|target(k)> = value(k) of
  !> configuration i, k = first(i) .. first(i + 1) - 1.
  type :: T_Hops
    integer(I32), allocatable :: first(:)
    integer(I32), allocatable :: target(:)
    complex(R64), allocatable :: value(:)
    complex(R64), allocatable :: diag(:)
  end type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Coeffs_Hubbard_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Coeffs
    use M_Coeffs_Hubbard_NoSpinSym
    use M_Coeffs_Hubbard_PlusSpinSym
    use M_Coeffs_Hubbard_MinusSpinSym

    call Say_Fabricate("coeffs.hubbard")

    !------------------------------------
    ! many-body setting
    !------------------------------------

    if (Method_Mb_nBodyTypes .ne. 2) error stop "coeffs.hubbard needs two body types (spin up and spin down)"
    if (any(Method_Mb_bodyStatistics .ne. 'f')) error stop "coeffs.hubbard needs fermions (bosons: coeffs.generic)"
    if (Method_Mb_OrbBased_nOrbs(1) .ne. Method_Mb_OrbBased_nOrbs(2)) &
      error stop "coeffs.hubbard needs the same number of orbitals (the lattice sites) for both spins"

    Coeffs_Hubbard_nSites = Method_Mb_OrbBased_nOrbs(1)
    if (Coeffs_Hubbard_nSites > 62) error stop "coeffs.hubbard supports at most 62 sites"

    call FillBinomials(binomials, Coeffs_Hubbard_nSites)
    Coeffs_Hubbard_nCoeffsUP = NumberOfConfigurations(Method_Mb_nBodies(1))
    Coeffs_Hubbard_nCoeffsDN = NumberOfConfigurations(Method_Mb_nBodies(2))

    !------------------------------------
    ! set procedure pointers
    !------------------------------------

    Coeffs_Setup => Setup
    Coeffs_ApplyH1FillRdm1 => ApplyH1FillRdm1
    Coeffs_ApplyH2FillRdm2 => ApplyH2FillRdm2
    Coeffs_ApplyExcitation => ApplyExcitation
    Coeffs_ConfigurationsFromIndex => ConfigurationsFromIndex
    Coeffs_IndexFromConfigurations => IndexFromConfigurations

    !------------------------------------
    ! branch on the spin symmetry
    !------------------------------------

    if (Json_GetExistence("coeffs.hubbard.noSpinSym")) then
      call Coeffs_Hubbard_NoSpinSym_Fabricate

    else if (Json_GetExistence("coeffs.hubbard.plusSpinSym")) then
      call Coeffs_Hubbard_PlusSpinSym_Fabricate

    else if (Json_GetExistence("coeffs.hubbard.minusSpinSym")) then
      call Coeffs_Hubbard_MinusSpinSym_Fabricate

    else
      error stop "coeffs.hubbard is missing one of: noSpinSym, plusSpinSym, minusSpinSym"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Configurations of both spins and their single excitations. Repeated calls
!> rebuild the tables.
  subroutine Setup
    use M_Utils_Say
    use M_Method_Mb

    call Say_Setup("coeffs.hubbard")

    call FillConfigurations(Coeffs_Hubbard_bitcodesUP, Method_Mb_nBodies(1))
    call FillConfigurations(Coeffs_Hubbard_bitcodesDN, Method_Mb_nBodies(2))

    call FillSingles(Coeffs_Hubbard_singlesUP, Coeffs_Hubbard_bitcodesUP)
    call FillSingles(Coeffs_Hubbard_singlesDN, Coeffs_Hubbard_bitcodesDN)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> H1 = sum_ij h1(i,j) a+_i a_j added to dCoeffs_ (h1_ of the call), and the
!> 1-RDM. A sector needs a spin-independent h1.
  subroutine ApplyH1FillRdm1(coeffs, rdm1_, dCoeffs_, h1_)
    complex(R64), intent(in), contiguous, target                 :: coeffs(:)
    complex(R64), intent(out), allocatable, optional             :: rdm1_(:, :)
    complex(R64), intent(inout), contiguous, target, optional    :: dCoeffs_(:)
    complex(R64), intent(in), contiguous, optional               :: h1_(:, :)

    complex(R64), allocatable :: c(:), dc(:)

    if (present(dCoeffs_) .and. .not. present(h1_)) error stop "Coeffs_ApplyH1FillRdm1 (coeffs.hubbard): dCoeffs_ requires h1_"

    if (Coeffs_Hubbard_parity .eq. 0) then
      if (present(dCoeffs_)) call AddH1(dCoeffs_, coeffs, h1_)
      if (present(rdm1_)) call FillRdm1(rdm1_, coeffs)
      return
    end if

    call ExpandSector(c, coeffs)

    if (present(dCoeffs_)) then
      call CheckSpinIndependent(h1_)
      allocate (dc(size(c)))
      dc(:) = 0.0_R64
      call AddH1(dc, c, h1_)
      call AddPackedSector(dCoeffs_, dc)
    end if

    if (present(rdm1_)) call FillRdm1(rdm1_, c)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The on-site interaction of h2_ added to dCoeffs_ (stops if h2_ is not an
!> on-site interaction), and the 2-RDM.
  subroutine ApplyH2FillRdm2(coeffs, rdm2_, dCoeffs_, h2_)
    complex(R64), intent(in), contiguous, target                 :: coeffs(:)
    complex(R64), intent(out), allocatable, optional             :: rdm2_(:, :, :, :)
    complex(R64), intent(inout), contiguous, target, optional    :: dCoeffs_(:)
    complex(R64), intent(in), contiguous, optional               :: h2_(:, :, :, :)

    complex(R64), allocatable :: c(:), dc(:)

    if (present(dCoeffs_) .and. .not. present(h2_)) error stop "Coeffs_ApplyH2FillRdm2 (coeffs.hubbard): dCoeffs_ requires h2_"

    if (Coeffs_Hubbard_parity .eq. 0) then
      if (present(dCoeffs_)) call AddH2(dCoeffs_, coeffs, h2_)
      if (present(rdm2_)) call FillRdm2(rdm2_, coeffs)
      return
    end if

    call ExpandSector(c, coeffs)

    if (present(dCoeffs_)) then
      allocate (dc(size(c)))
      dc(:) = 0.0_R64
      call AddH2(dc, c, h2_)
      call AddPackedSector(dCoeffs_, dc)
    end if

    if (present(rdm2_)) call FillRdm2(rdm2_, c)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> a+_creates(1) .. a+_creates(n) a_destroys(n) .. a_destroys(1) on the spin
!> bt (orbital indices of that spin), in place. Product space only: an
!> excitation of one spin leaves the spin-exchange sectors.
  subroutine ApplyExcitation(coeffs, creates, destroys, bt)
    complex(R64), intent(inout), contiguous :: coeffs(:)
    integer(I32), intent(in), contiguous :: creates(:)
    integer(I32), intent(in), contiguous :: destroys(:)
    integer(I32), intent(in) :: bt

    complex(R64), allocatable :: source(:)
    integer(I32) :: nUP, nDN, i, j, k, sgn
    integer(I64) :: y

    if (Coeffs_Hubbard_parity .ne. 0) &
      error stop "coeffs.hubbard: an excitation of one spin leaves the spin-exchange sector (use noSpinSym)"
    if (size(creates) .ne. size(destroys)) error stop "coeffs.hubbard: an excitation must conserve the particle number"
    if (any(creates < 1 .or. creates > Coeffs_Hubbard_nSites) .or. any(destroys < 1 .or. destroys > Coeffs_Hubbard_nSites)) &
      error stop "coeffs.hubbard: orbital index of an excitation out of range"

    nUP = Coeffs_Hubbard_nCoeffsUP
    nDN = Coeffs_Hubbard_nCoeffsDN

    allocate (source, source=coeffs)
    coeffs(:) = 0.0_R64

    if (bt .eq. 1) then
      do i = 1, nUP
        call ExciteBits(y, sgn, creates, destroys, Coeffs_Hubbard_bitcodesUP(i))
        if (sgn .eq. 0) cycle
        j = ConfigurationRank(y)
        do k = 0, nDN - 1
          coeffs(j + k * nUP) = real(sgn, kind=R64) * source(i + k * nUP)
        end do
      end do

    else if (bt .eq. 2) then
      do i = 1, nDN
        call ExciteBits(y, sgn, creates, destroys, Coeffs_Hubbard_bitcodesDN(i))
        if (sgn .eq. 0) cycle
        j = ConfigurationRank(y)
        coeffs((j - 1) * nUP + 1:j * nUP) = real(sgn, kind=R64) * source((i - 1) * nUP + 1:i * nUP)
      end do

    else
      error stop "coeffs.hubbard: the body type of an excitation must be 1 (up) or 2 (down)"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Configurations of the stored coefficient iCoeff: (iUp, iDn) in the product
!> space, (lo, hi) of the sector basis state in a sector.
  pure subroutine ConfigurationsFromIndex(configurations, iCoeff)
    integer(I32), intent(out), contiguous :: configurations(:)
    integer(I32), intent(in) :: iCoeff

    integer(I32) :: hi

    if (Coeffs_Hubbard_parity .eq. 0) then
      configurations(1) = mod(iCoeff - 1, Coeffs_Hubbard_nCoeffsUP) + 1
      configurations(2) = (iCoeff - 1) / Coeffs_Hubbard_nCoeffsUP + 1
      return
    end if

    ! column hi of the upper triangle ends at lo = hi (plus) or hi - 1 (minus)
    if (Coeffs_Hubbard_parity .eq. 1) then
      hi = 1
      do while (PackedIndex(hi, hi) < iCoeff)
        hi = hi + 1
      end do
    else
      hi = 2
      do while (PackedIndex(hi - 1, hi) < iCoeff)
        hi = hi + 1
      end do
    end if
    configurations(2) = hi
    configurations(1) = iCoeff - PackedIndex(1, hi) + 1

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Inverse of ConfigurationsFromIndex; in a sector both orders of the pair
!> give the same coefficient, and 0 for the pair (i, i) of the minus sector,
!> which has no coefficient.
  pure subroutine IndexFromConfigurations(iCoeff, configurations)
    integer(I32), intent(out) :: iCoeff
    integer(I32), intent(in), contiguous :: configurations(:)

    integer(I32) :: lo, hi

    if (Coeffs_Hubbard_parity .eq. 0) then
      iCoeff = configurations(1) + (configurations(2) - 1) * Coeffs_Hubbard_nCoeffsUP
      return
    end if

    lo = minval(configurations(1:2))
    hi = maxval(configurations(1:2))
    iCoeff = 0
    if (lo .ne. hi .or. Coeffs_Hubbard_parity .eq. 1) iCoeff = PackedIndex(lo, hi)

  end subroutine

!=============================================================================
! Hamiltonian application and RDMs on the product space c(iUp + (iDn - 1) nUP)
!=============================================================================

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> dc = dc + H1 c. The entry a+_c a_d |i> = factor |t> of configuration i
!> gives <i| a+_d a_c |t> = factor, so (H1 c)(i) = sum over the entries of i
!> of h1(d, c) factor c(t).
  subroutine AddH1(dc, c, h1)
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous :: dc(:)
    complex(R64), intent(in), contiguous :: c(:)
    complex(R64), intent(in), contiguous :: h1(:, :)

    type(T_Hops) :: hops
    integer(I32) :: n, oUp, oDn

    n = Coeffs_Hubbard_nSites
    oUp = Method_Mb_OrbBased_nOrbsStart(1) - 1
    oDn = Method_Mb_OrbBased_nOrbsStart(2) - 1

    call FillHops(hops, Coeffs_Hubbard_singlesUP, h1(oUp + 1:oUp + n, oUp + 1:oUp + n))
    call AddH1Up(dc, c, hops)

    call FillHops(hops, Coeffs_Hubbard_singlesDN, h1(oDn + 1:oDn + n, oDn + 1:oDn + n))
    call AddH1Down(dc, c, hops)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The spin-up part of H1: within every column iDn of c(iUp, iDn).
  subroutine AddH1Up(dc, c, up)
    complex(R64), intent(inout), contiguous :: dc(:)
    complex(R64), intent(in), contiguous :: c(:)
    type(T_Hops), intent(in) :: up

    complex(R64) :: acc
    integer(I32) :: nUP, iUp, iDn, k, base

    nUP = Coeffs_Hubbard_nCoeffsUP

    !$omp parallel do default(shared) private(iDn, iUp, k, acc, base) schedule(static)
    do iDn = 1, Coeffs_Hubbard_nCoeffsDN
      base = (iDn - 1) * nUP
      do iUp = 1, nUP
        acc = up % diag(iUp) * c(base + iUp)
        do k = up % first(iUp), up % first(iUp + 1) - 1
          acc = acc + up % value(k) * c(base + up % target(k))
        end do
        dc(base + iUp) = dc(base + iUp) + acc
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The spin-down part of H1: whole columns of c(iUp, iDn).
  subroutine AddH1Down(dc, c, dn)
    complex(R64), intent(inout), contiguous :: dc(:)
    complex(R64), intent(in), contiguous :: c(:)
    type(T_Hops), intent(in) :: dn

    integer(I32) :: nUP, iDn, k, base, t

    nUP = Coeffs_Hubbard_nCoeffsUP

    !$omp parallel do default(shared) private(iDn, k, base, t) schedule(static)
    do iDn = 1, Coeffs_Hubbard_nCoeffsDN
      base = (iDn - 1) * nUP
      dc(base + 1:base + nUP) = dc(base + 1:base + nUP) + dn % diag(iDn) * c(base + 1:base + nUP)
      do k = dn % first(iDn), dn % first(iDn + 1) - 1
        t = (dn % target(k) - 1) * nUP
        dc(base + 1:base + nUP) = dc(base + 1:base + nUP) + dn % value(k) * c(t + 1:t + nUP)
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The H1 of one spin (h1s the block of h1 of that spin) from its single
!> excitations: the entries with c = d sum to the diagonal, and of the others
!> only the nonzero elements are kept (the hops of the lattice).
  subroutine FillHops(hops, singles, h1s)
    type(T_Hops), intent(out) :: hops
    type(T_Coeffs_Hubbard_Singles), intent(in) :: singles
    complex(R64), intent(in), contiguous :: h1s(:, :)

    complex(R64) :: v
    integer(I32) :: nConf, i, k, m

    nConf = size(singles % first) - 1
    allocate (hops % first(nConf + 1), hops % diag(nConf))
    allocate (hops % target(size(singles % target)), hops % value(size(singles % target)))

    m = 0
    do i = 1, nConf
      hops % first(i) = m + 1
      hops % diag(i) = 0.0_R64
      do k = singles % first(i), singles % first(i + 1) - 1
        v = singles % factor(k) * h1s(singles % destroy(k), singles % create(k))
        if (singles % create(k) .eq. singles % destroy(k)) then
          hops % diag(i) = hops % diag(i) + v
        else if (abs(v) > 0.0_R64) then
          m = m + 1
          hops % target(m) = singles % target(k)
          hops % value(m) = v
        end if
      end do
    end do
    hops % first(nConf + 1) = m + 1

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> dc = dc + H2 c with H2 = sum_s U_s n_{s,up} n_{s,dn} (diagonal).
  subroutine AddH2(dc, c, h2)
    complex(R64), intent(inout), contiguous :: dc(:)
    complex(R64), intent(in), contiguous :: c(:)
    complex(R64), intent(in), contiguous :: h2(:, :, :, :)

    complex(R64), allocatable :: u(:)
    complex(R64) :: w
    integer(I64) :: both
    integer(I32) :: nUP, nDN, iUp, iDn, s, idx

    call FillOnSiteInteraction(u, h2)

    nUP = Coeffs_Hubbard_nCoeffsUP
    nDN = Coeffs_Hubbard_nCoeffsDN

    !$omp parallel do default(shared) private(iDn, iUp, both, w, s, idx) schedule(static)
    do iDn = 1, nDN
      do iUp = 1, nUP
        ! doubly occupied sites
        both = iand(Coeffs_Hubbard_bitcodesUP(iUp), Coeffs_Hubbard_bitcodesDN(iDn))
        if (both .eq. 0) cycle
        w = (0.0_R64, 0.0_R64)
        do while (both .ne. 0)
          s = trailz(both)
          w = w + u(s + 1)
          both = ibclr(both, s)
        end do
        idx = iUp + (iDn - 1) * nUP
        dc(idx) = dc(idx) + w * c(idx)
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> U_s of the on-site interaction that h2 represents. For up orbitals a, b
!> and down orbitals g, d coeffs.generic couples the excitations a+_a a_b
!> (up) and a+_g a_d (down) by
!>   v = h2(b,d,a,g) - h2(d,b,a,g) - h2(b,d,g,a) + h2(d,b,g,a)
!> (spin-orbital indices); U_s = v at a = b = g = d = s. Stops unless v
!> vanishes for every other (a, b, g, d) and the antisymmetrized interaction
!> within each spin vanishes (below 1e-10 times the largest of these
!> combinations, or 1): then this backend and coeffs.generic apply the same
!> H2. Only these combinations are read, O(n^4) per call.
  subroutine FillOnSiteInteraction(u, h2)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: u(:)
    complex(R64), intent(in), contiguous :: h2(:, :, :, :)

    integer(I32) :: n, oUp, oDn, a, b, g, d, bt, o
    integer(I32) :: i1, i2, j1, j2
    complex(R64) :: v
    real(R64) :: tol, offSite, sameSpin

    n = Coeffs_Hubbard_nSites
    oUp = Method_Mb_OrbBased_nOrbsStart(1) - 1
    oDn = Method_Mb_OrbBased_nOrbsStart(2) - 1

    allocate (u(n))
    offSite = 0.0_R64

    do d = 1, n
      do g = 1, n
        do b = 1, n
          do a = 1, n
            v = h2(oUp + b, oDn + d, oUp + a, oDn + g) - h2(oDn + d, oUp + b, oUp + a, oDn + g) &
                - h2(oUp + b, oDn + d, oDn + g, oUp + a) + h2(oDn + d, oUp + b, oDn + g, oUp + a)
            if (a .eq. b .and. a .eq. g .and. a .eq. d) then
              u(a) = v
            else
              offSite = max(offSite, abs(v))
            end if
          end do
        end do
      end do
    end do

    sameSpin = 0.0_R64
    do bt = 1, 2
      o = Method_Mb_OrbBased_nOrbsStart(bt) - 1
      do j2 = o + 2, o + n
        do j1 = o + 1, j2 - 1
          do i2 = o + 2, o + n
            do i1 = o + 1, i2 - 1
              v = h2(j1, j2, i1, i2) - h2(j2, j1, i1, i2) - h2(j1, j2, i2, i1) + h2(j2, j1, i2, i1)
              sameSpin = max(sameSpin, abs(v))
            end do
          end do
        end do
      end do
    end do

    tol = 1e-10_R64 * max(1.0_R64, maxval(abs(u)), offSite, sameSpin)
    if (offSite > tol .or. sameSpin > tol) then
      print *, "coeffs.hubbard: largest off-site interaction ", offSite, ", largest same-spin interaction ", sameSpin
      error stop "coeffs.hubbard: h2 is not an on-site interaction of the orbitals (lattice sites, static orbitals: tdci)"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> rdm1(p, q) = <a+_q a_p>: the entry a+_c a_d |i> = factor |t> adds
!> factor <t|i> (overlap over the other spin) to rdm1(d, c).
  subroutine FillRdm1(rdm1, c)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm1(:, :)
    complex(R64), intent(in), contiguous :: c(:)

    complex(R64), allocatable :: block(:, :)
    integer(I32) :: n, oUp, oDn

    n = Coeffs_Hubbard_nSites
    oUp = Method_Mb_OrbBased_nOrbsStart(1) - 1
    oDn = Method_Mb_OrbBased_nOrbsStart(2) - 1

    allocate (rdm1(Method_Mb_OrbBased_nOrbsSum, Method_Mb_OrbBased_nOrbsSum))
    rdm1(:, :) = 0.0_R64

    call FillRdm1Block(block, c, 1, Coeffs_Hubbard_singlesUP)
    rdm1(oUp + 1:oUp + n, oUp + 1:oUp + n) = block
    call FillRdm1Block(block, c, 2, Coeffs_Hubbard_singlesDN)
    rdm1(oDn + 1:oDn + n, oDn + 1:oDn + n) = block

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The 1-RDM block of spin bt (orbital indices of that spin) from its single
!> excitations.
  subroutine FillRdm1Block(block, c, bt, singles)
    complex(R64), intent(out), allocatable :: block(:, :)
    complex(R64), intent(in), contiguous :: c(:)
    integer(I32), intent(in) :: bt
    type(T_Coeffs_Hubbard_Singles), intent(in) :: singles

    complex(R64), allocatable :: strings(:, :), blockThread(:, :)
    integer(I32) :: n, i, k, cr, de

    n = Coeffs_Hubbard_nSites
    call FillStrings(strings, c, bt)

    allocate (block(n, n))
    block(:, :) = 0.0_R64

    !$omp parallel default(shared) private(i, k, cr, de, blockThread)
    allocate (blockThread(n, n))
    blockThread(:, :) = 0.0_R64
    !$omp do schedule(static)
    do i = 1, size(strings, 2)
      do k = singles % first(i), singles % first(i + 1) - 1
        cr = singles % create(k)
        de = singles % destroy(k)
        blockThread(de, cr) = blockThread(de, cr) &
                              + singles % factor(k) * dot_product(strings(:, singles % target(k)), strings(:, i))
      end do
    end do
    !$omp end do
    !$omp critical
    block = block + blockThread
    !$omp end critical
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> rdm2(p,q,r,s) = <a+_r a+_s a_q a_p>. Within one spin
!>   <a+_r a+_s a_q a_p> = <a+_r a_p a+_s a_q> - delta(p,s) <a+_r a_q>
!> from chains of two single excitations; for p, r up and q, s down
!>   X(p,q,r,s) = <(a+_r a_p)_up (a+_s a_q)_dn>
!> gives rdm2(p,q,r,s) = rdm2(q,p,s,r) = X and, antisymmetric in the
!> exchange of an up and a down index (the convention of coeffs.generic),
!> rdm2(p,q,s,r) = rdm2(q,p,r,s) = -X.
  subroutine FillRdm2(rdm2, c)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous :: c(:)

    complex(R64), allocatable :: rdm1(:, :), x(:, :, :, :)
    integer(I32) :: n, bt, o, oUp, oDn, p, q, r, s

    n = Coeffs_Hubbard_nSites
    oUp = Method_Mb_OrbBased_nOrbsStart(1) - 1
    oDn = Method_Mb_OrbBased_nOrbsStart(2) - 1

    allocate (rdm2(Method_Mb_OrbBased_nOrbsSum, Method_Mb_OrbBased_nOrbsSum, &
                   Method_Mb_OrbBased_nOrbsSum, Method_Mb_OrbBased_nOrbsSum))
    rdm2(:, :, :, :) = 0.0_R64

    call FillRdm1(rdm1, c)

    ! within one spin
    do bt = 1, 2
      o = Method_Mb_OrbBased_nOrbsStart(bt) - 1
      if (bt .eq. 1) call FillSameSpinProducts(x, c, bt, Coeffs_Hubbard_singlesUP)
      if (bt .eq. 2) call FillSameSpinProducts(x, c, bt, Coeffs_Hubbard_singlesDN)
      do s = 1, n
        do r = 1, n
          do q = 1, n
            do p = 1, n
              rdm2(o + p, o + q, o + r, o + s) = x(p, q, r, s)
              if (p .eq. s) rdm2(o + p, o + q, o + r, o + s) = rdm2(o + p, o + q, o + r, o + s) - rdm1(o + q, o + r)
            end do
          end do
        end do
      end do
    end do

    ! between the spins
    call FillUpDownProducts(x, c, Coeffs_Hubbard_singlesUP, Coeffs_Hubbard_singlesDN)
    do s = 1, n
      do r = 1, n
        do q = 1, n
          do p = 1, n
            rdm2(oUp + p, oDn + q, oUp + r, oDn + s) = x(p, q, r, s)
            rdm2(oDn + q, oUp + p, oDn + s, oUp + r) = x(p, q, r, s)
            rdm2(oUp + p, oDn + q, oDn + s, oUp + r) = -x(p, q, r, s)
            rdm2(oDn + q, oUp + p, oUp + r, oDn + s) = -x(p, q, r, s)
          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> x(p,q,r,s) = <a+_r a_p a+_s a_q> within spin bt: the entry a+_s a_q |i> =
!> f1 |m> followed by the entry a+_r a_p |m> = f2 |t> adds f1 f2 <t|i>.
  subroutine FillSameSpinProducts(x, c, bt, singles)
    complex(R64), intent(out), allocatable :: x(:, :, :, :)
    complex(R64), intent(in), contiguous :: c(:)
    integer(I32), intent(in) :: bt
    type(T_Coeffs_Hubbard_Singles), intent(in) :: singles

    complex(R64), allocatable :: strings(:, :), xThread(:, :, :, :)
    integer(I32) :: n, i, k1, k2, m, p, q, r, s

    n = Coeffs_Hubbard_nSites
    call FillStrings(strings, c, bt)

    allocate (x(n, n, n, n))
    x(:, :, :, :) = 0.0_R64

    !$omp parallel default(shared) private(i, k1, k2, m, p, q, r, s, xThread)
    allocate (xThread(n, n, n, n))
    xThread(:, :, :, :) = 0.0_R64
    !$omp do schedule(dynamic)
    do i = 1, size(strings, 2)
      do k1 = singles % first(i), singles % first(i + 1) - 1
        m = singles % target(k1)
        q = singles % destroy(k1)
        s = singles % create(k1)
        do k2 = singles % first(m), singles % first(m + 1) - 1
          p = singles % destroy(k2)
          r = singles % create(k2)
          xThread(p, q, r, s) = xThread(p, q, r, s) + singles % factor(k1) * singles % factor(k2) &
                                * dot_product(strings(:, singles % target(k2)), strings(:, i))
        end do
      end do
    end do
    !$omp end do
    !$omp critical
    x = x + xThread
    !$omp end critical
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> x(p,q,r,s) = <(a+_r a_p)_up (a+_s a_q)_dn>: the entries a+_s a_q |iDn> =
!> fd |mDn> and a+_r a_p |iUp> = fu |tUp> add fu fd conj(c(tUp, mDn))
!> c(iUp, iDn).
  subroutine FillUpDownProducts(x, c, up, dn)
    complex(R64), intent(out), allocatable :: x(:, :, :, :)
    complex(R64), intent(in), contiguous :: c(:)
    type(T_Coeffs_Hubbard_Singles), intent(in) :: up
    type(T_Coeffs_Hubbard_Singles), intent(in) :: dn

    complex(R64), allocatable :: xThread(:, :, :, :)
    integer(I32) :: n, nUP, iUp, iDn, ku, kd, mBase, iBase, p, q, r, s

    n = Coeffs_Hubbard_nSites
    nUP = Coeffs_Hubbard_nCoeffsUP

    allocate (x(n, n, n, n))
    x(:, :, :, :) = 0.0_R64

    !$omp parallel default(shared) private(iUp, iDn, ku, kd, mBase, iBase, p, q, r, s, xThread)
    allocate (xThread(n, n, n, n))
    xThread(:, :, :, :) = 0.0_R64
    !$omp do schedule(dynamic)
    do iDn = 1, Coeffs_Hubbard_nCoeffsDN
      iBase = (iDn - 1) * nUP
      do kd = dn % first(iDn), dn % first(iDn + 1) - 1
        mBase = (dn % target(kd) - 1) * nUP
        q = dn % destroy(kd)
        s = dn % create(kd)
        do iUp = 1, nUP
          if (.not. (abs(c(iBase + iUp)) > 0.0_R64)) cycle
          do ku = up % first(iUp), up % first(iUp + 1) - 1
            p = up % destroy(ku)
            r = up % create(ku)
            xThread(p, q, r, s) = xThread(p, q, r, s) + up % factor(ku) * dn % factor(kd) &
                                  * conjg(c(mBase + up % target(ku))) * c(iBase + iUp)
          end do
        end do
      end do
    end do
    !$omp end do
    !$omp critical
    x = x + xThread
    !$omp end critical
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> strings(:, i) = the coefficients of configuration i of spin bt over the
!> configurations of the other spin (a row of c(iUp, iDn) for up, a column
!> for down), contiguous for the overlaps of the RDMs.
  subroutine FillStrings(strings, c, bt)
    complex(R64), intent(out), allocatable :: strings(:, :)
    complex(R64), intent(in), contiguous :: c(:)
    integer(I32), intent(in) :: bt

    if (bt .eq. 1) then
      strings = transpose(reshape(c, [Coeffs_Hubbard_nCoeffsUP, Coeffs_Hubbard_nCoeffsDN]))
    else
      strings = reshape(c, [Coeffs_Hubbard_nCoeffsUP, Coeffs_Hubbard_nCoeffsDN])
    end if

  end subroutine

!=============================================================================
! spin-exchange sectors
!=============================================================================

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Product-space vector c(iUp + (iDn - 1) n) of the sector vector packed:
!> c(lo, hi) = packed / sqrt(2), c(hi, lo) = parity packed / sqrt(2) for
!> lo < hi, c(i, i) = packed (plus sector) or 0 (minus sector).
  subroutine ExpandSector(c, packed)
    complex(R64), intent(out), allocatable :: c(:)
    complex(R64), intent(in), contiguous :: packed(:)

    integer(I32) :: m, lo, hi
    complex(R64) :: v

    m = Coeffs_Hubbard_nCoeffsUP
    allocate (c(m * m))
    c(:) = 0.0_R64

    ! iteration hi writes column hi and row hi of the upper part only
    !$omp parallel do default(shared) private(hi, lo, v) schedule(dynamic, 16)
    do hi = 1, m
      do lo = 1, hi - 1
        v = packed(PackedIndex(lo, hi)) / sqrt(2.0_R64)
        c(lo + (hi - 1) * m) = v
        c(hi + (lo - 1) * m) = real(Coeffs_Hubbard_parity, kind=R64) * v
      end do
      if (Coeffs_Hubbard_parity .eq. 1) c(hi + (hi - 1) * m) = packed(PackedIndex(hi, hi))
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> packed = packed + the components of the product-space vector dc on the
!> orthonormal sector basis, (dc(lo, hi) + parity dc(hi, lo)) / sqrt(2) and
!> dc(i, i) (plus sector).
  subroutine AddPackedSector(packed, dc)
    complex(R64), intent(inout), contiguous :: packed(:)
    complex(R64), intent(in), contiguous :: dc(:)

    integer(I32) :: m, lo, hi, idx

    m = Coeffs_Hubbard_nCoeffsUP

    !$omp parallel do default(shared) private(hi, lo, idx) schedule(dynamic, 16)
    do hi = 1, m
      do lo = 1, hi - 1
        idx = PackedIndex(lo, hi)
        packed(idx) = packed(idx) + (dc(lo + (hi - 1) * m) + real(Coeffs_Hubbard_parity, kind=R64) * dc(hi + (lo - 1) * m)) &
                      / sqrt(2.0_R64)
      end do
      if (Coeffs_Hubbard_parity .eq. 1) then
        idx = PackedIndex(hi, hi)
        packed(idx) = packed(idx) + dc(hi + (hi - 1) * m)
      end if
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Stored index of the sector basis state of the configurations lo <= hi
!> (plus sector) or lo < hi (minus sector), columns of the upper triangle.
  pure function PackedIndex(lo, hi) result(idx)
    integer(I32), intent(in) :: lo
    integer(I32), intent(in) :: hi
    integer(I32) :: idx

    if (Coeffs_Hubbard_parity .eq. 1) then
      idx = int((int(hi - 1, I64) * hi) / 2, I32) + lo
    else
      idx = int((int(hi - 2, I64) * (hi - 1)) / 2, I32) + lo
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The sectors are invariant only under an h1 that is the same for both spins.
  subroutine CheckSpinIndependent(h1)
    use M_Method_Mb_OrbBased

    complex(R64), intent(in), contiguous :: h1(:, :)

    integer(I32) :: n, oUp, oDn
    real(R64) :: difference

    n = Coeffs_Hubbard_nSites
    oUp = Method_Mb_OrbBased_nOrbsStart(1) - 1
    oDn = Method_Mb_OrbBased_nOrbsStart(2) - 1

    difference = maxval(abs(h1(oUp + 1:oUp + n, oUp + 1:oUp + n) - h1(oDn + 1:oDn + n, oDn + 1:oDn + n)))
    if (difference > 1e-12_R64 * max(1.0_R64, maxval(abs(h1)))) then
      print *, "coeffs.hubbard: largest difference of h1 between the spins ", difference
      error stop "coeffs.hubbard: the spin-exchange sectors need an h1 that is the same for both spins (use noSpinSym)"
    end if

  end subroutine

!=============================================================================
! configurations and bit operations
!=============================================================================

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillBinomials(b, n)
    integer(I64), intent(out), allocatable :: b(:, :)
    integer(I32), intent(in) :: n

    integer(I32) :: i, k

    allocate (b(0:n, 0:n))
    b(:, :) = 0
    do i = 0, n
      b(i, 0) = 1
      do k = 1, i
        b(i, k) = b(i - 1, k - 1) + b(i - 1, k)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function NumberOfConfigurations(nParticles) result(nConfigurations)
    integer(I32), intent(in) :: nParticles
    integer(I32) :: nConfigurations

    if (nParticles < 0 .or. nParticles > Coeffs_Hubbard_nSites) &
      error stop "coeffs.hubbard: more particles of one spin than orbitals"
    if (binomials(Coeffs_Hubbard_nSites, nParticles) > huge(nConfigurations)) &
      error stop "coeffs.hubbard: more than 2^31 - 1 configurations of one spin"

    nConfigurations = int(binomials(Coeffs_Hubbard_nSites, nParticles), I32)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> All bit patterns of nParticles bits below bit nSites, in increasing order:
!> from the smallest, the next larger pattern with as many bits (Gosper). Their
!> colex ranks must count 1, 2, ...
  subroutine FillConfigurations(bitcodes, nParticles)
    integer(I64), intent(out), allocatable :: bitcodes(:)
    integer(I32), intent(in) :: nParticles

    integer(I32) :: iConf
    integer(I64) :: x, lowest, ripple

    allocate (bitcodes(NumberOfConfigurations(nParticles)))

    x = maskr(nParticles, I64)
    do iConf = 1, size(bitcodes)
      bitcodes(iConf) = x
      if (ConfigurationRank(x) .ne. iConf) error stop "coeffs.hubbard: the rank of a configuration differs from its position"
      if (iConf .eq. size(bitcodes)) exit
      lowest = iand(x, -x)
      ripple = x + lowest
      x = ior(ishft(ieor(ripple, x), -2) / lowest, ripple)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Single excitations a+_c a_d (d occupied, c empty or c = d) of every
!> configuration: N (nSites - N) + N entries per configuration.
  subroutine FillSingles(singles, bitcodes)
    type(T_Coeffs_Hubbard_Singles), intent(out) :: singles
    integer(I64), intent(in), contiguous :: bitcodes(:)

    integer(I32) :: nConf, nParticles, nPer, iConf, k, c, d
    integer(I64) :: x, nEntries

    nConf = size(bitcodes)
    nParticles = popcnt(bitcodes(1))
    nPer = nParticles * (Coeffs_Hubbard_nSites - nParticles) + nParticles
    nEntries = int(nConf, I64) * nPer
    if (nEntries > huge(nConf)) error stop "coeffs.hubbard: more than 2^31 - 1 single excitations of one spin"

    allocate (singles % first(nConf + 1))
    allocate (singles % target(nEntries), singles % create(nEntries), singles % destroy(nEntries), &
              singles % factor(nEntries))

    !$omp parallel do default(shared) private(iConf, k, c, d, x) schedule(static)
    do iConf = 1, nConf
      singles % first(iConf) = (iConf - 1) * nPer + 1
      k = (iConf - 1) * nPer
      x = bitcodes(iConf)
      do d = 1, Coeffs_Hubbard_nSites
        if (.not. btest(x, d - 1)) cycle
        do c = 1, Coeffs_Hubbard_nSites
          if (c .ne. d .and. btest(x, c - 1)) cycle
          k = k + 1
          singles % target(k) = ConfigurationRank(ibset(ibclr(x, d - 1), c - 1))
          singles % create(k) = c
          singles % destroy(k) = d
          singles % factor(k) = HopSign(x, c, d)
        end do
      end do
    end do
    !$omp end parallel do
    singles % first(nConf + 1) = nConf * nPer + 1

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Index of the configuration x among the configurations of its number of
!> particles: one plus its colex rank sum_t binom(position of bit t, t).
  pure function ConfigurationRank(x) result(iConf)
    integer(I64), intent(in) :: x
    integer(I32) :: iConf

    integer(I64) :: rest, r
    integer(I32) :: t, pos

    rest = x
    r = 0
    t = 0
    do while (rest .ne. 0)
      pos = trailz(rest)
      t = t + 1
      r = r + binomials(pos, t)
      rest = ibclr(rest, pos)
    end do
    iConf = int(r, I32) + 1

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Sign of a+_c a_d on x (d occupied): (-1)^(occupied orbitals strictly
!> between c and d).
  pure function HopSign(x, c, d) result(sgn)
    integer(I64), intent(in) :: x
    integer(I32), intent(in) :: c
    integer(I32), intent(in) :: d
    real(R64) :: sgn

    integer(I32) :: lo, hi

    lo = min(c, d)
    hi = max(c, d)
    sgn = 1.0_R64
    if (hi - lo .le. 1) return
    ! orbitals lo + 1 .. hi - 1 are the bits lo .. hi - 2
    if (mod(popcnt(ibits(x, lo, hi - lo - 1)), 2) .eq. 1) sgn = -1.0_R64

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> y = a+_creates(1) .. a+_creates(n) a_destroys(n) .. a_destroys(1) x up to
!> the sign sgn (0 if the excitation vanishes); every operator on orbital o
!> picks up (-1)^(occupied orbitals below o).
  pure subroutine ExciteBits(y, sgn, creates, destroys, x)
    integer(I64), intent(out) :: y
    integer(I32), intent(out) :: sgn
    integer(I32), intent(in), contiguous :: creates(:)
    integer(I32), intent(in), contiguous :: destroys(:)
    integer(I64), intent(in) :: x

    integer(I32) :: i, o

    y = x
    sgn = 1

    do i = 1, size(destroys)
      o = destroys(i)
      if (.not. btest(y, o - 1)) then
        sgn = 0
        return
      end if
      if (mod(popcnt(ibits(y, 0, o - 1)), 2) .eq. 1) sgn = -sgn
      y = ibclr(y, o - 1)
    end do

    do i = size(creates), 1, -1
      o = creates(i)
      if (btest(y, o - 1)) then
        sgn = 0
        return
      end if
      if (mod(popcnt(ibits(y, 0, o - 1)), 2) .eq. 1) sgn = -sgn
      y = ibset(y, o - 1)
    end do

  end subroutine

end submodule
