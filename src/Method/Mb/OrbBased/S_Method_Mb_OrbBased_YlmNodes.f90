! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased:S_Method_Mb_OrbBased) S_Method_Mb_OrbBased_YlmNodes
!-----------------------------------------------------------------------------
! Definite-m node path of the orbital-based methods (Ylm grids with
! method.mb.orbBased.mOfOrbs): the correlation operator and h2 on the
! Gauss-Legendre nodes in cos(theta), where all products are pointwise.
!
! Provides (declared in S_Method_Mb_OrbBased, which binds them):
!   - YlmNodesFastPathQ:          whether the node path applies
!   - SetupYlmNodes:              skeleton of pairs and target groups
!   - ApplyCorrelationOpYlmNodes: Method_Mb_OrbBased_ApplyCorrelationOp
!   - FillH2YlmNodes:             Method_Mb_OrbBased_FillH2
!
! A child of S_Method_Mb_OrbBased (nested submodule): FillCorrelationRdm2
! and DuplicateRestrictedSpinBlocks are reached by host association.
!-----------------------------------------------------------------------------

implicit none

!=============================================================================
! local data
!=============================================================================

!> Target group (m, body type) of the definite-m node path: its columns
!> (pair colPair, orbital colOrb), those of stored pairs first (nColStored,
!> the only ones entering h2), its orbitals, the offset of its stored
!> columns in the h2 columns of all groups and the 2-RDM weights
!> (columns x orbitals).
type :: T_TargetGroup
  integer(I32), allocatable :: colPair(:)
  integer(I32), allocatable :: colOrb(:)
  integer(I32), allocatable :: orbs(:)
  integer(I32) :: nColStored = 0
  integer(I32) :: colOffset = 0
  complex(R64), allocatable :: weights(:, :)
end type

!> Skeleton of the node path (SetupYlmNodes): ordered pairs (k2, l2), every
!> unordered pair followed by its swap, with M = m_k2 - m_l2 and the index of
!> the stored potential (negative for the swap: conjugate); the target groups
!> with their total stored (h2) and largest column count; the quadrature
!> weights.
integer(I32), allocatable :: nodePairK2(:), nodePairL2(:), nodePairM(:), nodePairPot(:)
integer(I32) :: nodeNPots = 0
type(T_TargetGroup), allocatable :: nodeGroups(:)
integer(I32) :: nodeNH2Cols = 0, nodeMaxCols = 0
real(R64), allocatable :: nodeWeights(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module function YlmNodesFastPathQ() result(res)
  !---------------------------------------------------------------------------
  ! Whether the definite-m node quadrature applies: Hermitian Ylm grid and a
  ! body-type-independent Ylm interaction within the node tables.
  !---------------------------------------------------------------------------
  use M_Utils_Json
  use M_Grid_Ylm
  use M_SysInteraction
  use M_SysInteraction_Ylm

  logical :: res

  res = .false.
  if (Grid_Ylm_cProductQ) return
  if (.not. Json_GetExistence("sysInteraction.ylm")) return
  if (.not. SysInteraction_bodyTypeIndependentQ) return
  if (SysInteraction_Ylm_lmax > Grid_Ylm_lmaxNodes) return
  if (SysInteraction_Ylm_mmax > Grid_Ylm_mmaxNodes) return

  res = .true.

end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module subroutine ApplyCorrelationOpYlmNodes(dOrbs, orbs, invRdm1, rdm2, time, h2_)
  !---------------------------------------------------------------------------
  ! ApplyCorrelationOp for orbitals of declared m on the Gauss-Legendre
  ! nodes in cos(theta), where all products are pointwise: the orbitals are
  ! transformed once, each pair potential once and dOrbs once at the end.
  ! The pair (l2,k2) carries M = m_k2 - m_l2; only the terms with
  ! m_k1 + M = m_j1 (correlation) and m_j1 + M = m_i1 (h2) are kept, the
  ! others vanish. Equal to ApplyCorrelationOp up to round-off (A_Method.md).
  !---------------------------------------------------------------------------
  use M_Grid
  use M_Grid_Ylm
  use M_Orbs

  complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
  complex(R64), intent(in), contiguous, target    :: orbs(:, :)
  complex(R64), intent(in), contiguous             :: invRdm1(:, :)
  complex(R64), intent(in), contiguous             :: rdm2(:, :, :, :)
  real(R64), intent(in)                  :: time
  complex(R64), intent(out), allocatable, optional :: h2_(:, :, :, :)

  complex(R64), allocatable :: rdm2spatial(:, :, :, :)
  complex(R64), allocatable, target :: orbNodes(:, :), dOrbNodes(:, :), potNodes(:, :)
  complex(R64), allocatable :: h2Cols(:, :), dOrb(:)
  complex(R64), contiguous, pointer :: nodes(:, :)

  integer(I32) :: nOS, nRad, nNodes, i1

  nOS = Orbs_nOrbsInState
  nRad = Grid_Ylm_nRadial
  nNodes = Grid_Ylm_nNodes

  call FillCorrelationRdm2(rdm2spatial, invRdm1, rdm2)
  call FillTargetGroupWeights(rdm2spatial)

  call FillOrbNodes(orbNodes, orbs)
  call FillPairPotentialsYlmNodes(potNodes, orbNodes, time)

  allocate (dOrbNodes, mold=orbNodes)
  dOrbNodes(:, :) = 0.0_R64

  if (present(h2_)) then
    call AccumulateTargetSumsYlmNodes(orbNodes, potNodes, h2Cols_=h2Cols, dOrbNodes_=dOrbNodes)
    call FillH2FromTargetColumns(h2_, h2Cols)
  else
    call AccumulateTargetSumsYlmNodes(orbNodes, potNodes, dOrbNodes_=dOrbNodes)
  end if

  ! Back to the Ylm representation
  !$omp parallel default(shared) private(i1, dOrb, nodes)
  allocate (dOrb(Grid_nPoints))
  !$omp do
  do i1 = 1, nOS
    nodes(1:nRad, 1:nNodes) => dOrbNodes(:, i1)
    call Grid_Ylm_FromNodes(dOrb, nodes, Method_Mb_OrbBased_mOfOrb(i1))
    dOrbs(:, i1) = dOrbs(:, i1) + dOrb(:)
  end do
  !$omp end do
  deallocate (dOrb)
  !$omp end parallel

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module subroutine FillH2YlmNodes(h2, orbs, time)
  !---------------------------------------------------------------------------
  ! FillH2 for orbitals of declared m by the node quadrature of
  ! ApplyCorrelationOpYlmNodes (pair potentials and h2 columns only):
  !   h2(i1,l2,k1,k2) = 1/2 <orb_i1| V_(l2,k2) |orb_k1>,  m_i1 = m_k1 + m_k2 - m_l2.
  !---------------------------------------------------------------------------
  complex(R64), intent(out), allocatable :: h2(:, :, :, :)
  complex(R64), intent(in), contiguous   :: orbs(:, :)
  real(R64), intent(in)                  :: time

  complex(R64), allocatable :: orbNodes(:, :), potNodes(:, :), h2Cols(:, :)

  call FillOrbNodes(orbNodes, orbs)
  call FillPairPotentialsYlmNodes(potNodes, orbNodes, time)
  call AccumulateTargetSumsYlmNodes(orbNodes, potNodes, h2Cols_=h2Cols)
  call FillH2FromTargetColumns(h2, h2Cols)

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module subroutine SetupYlmNodes
  !---------------------------------------------------------------------------
  ! Skeleton of the node path from the declared m of the stored orbitals: the
  ! pair list, the target groups (the blocks Orbs_blockOfOrb) of the
  ! selection rule m_k + M_iP = m_j with their columns (iP, k), and the
  ! quadrature weights on the node points.
  !---------------------------------------------------------------------------
  use M_Grid_Ylm
  use M_SysInteraction_Ylm
  use M_Orbs

  integer(I32), allocatable :: mOrb(:), btOrb(:), targetGroup(:, :), pairOfCol(:, :), orbOfCol(:, :)
  logical, allocatable :: storedQ(:, :)
  integer(I32) :: nOS, nPairs, iG, iP, j1, k1, k2, l2

  ! A repeated setup starts from an empty skeleton
  if (allocated(nodePairK2)) deallocate (nodePairK2, nodePairL2, nodePairM, nodePairPot)
  if (allocated(nodeGroups)) deallocate (nodeGroups)
  nodeNPots = 0
  nodeNH2Cols = 0
  nodeMaxCols = 0

  nOS = Orbs_nOrbsInState
  mOrb = Method_Mb_OrbBased_mOfOrb(1:nOS)
  btOrb = Method_Mb_OrbBased_bodyTypeOfOrb(1:nOS)

  ! Ordered pairs (k2, l2) of equal body type within the |m| cap of the
  ! potentials (pairs beyond it are left out, as in the generic operator)
  nPairs = count(spread(btOrb, 1, nOS) .eq. spread(btOrb, 2, nOS) .and. &
                 abs(spread(mOrb, 1, nOS) - spread(mOrb, 2, nOS)) <= SysInteraction_Ylm_mmax)
  allocate (nodePairK2(nPairs), nodePairL2(nPairs), nodePairM(nPairs), nodePairPot(nPairs))

  iP = 0
  do l2 = 1, nOS
    do k2 = l2, nOS
      if (btOrb(l2) .ne. btOrb(k2)) cycle
      if (abs(mOrb(k2) - mOrb(l2)) > SysInteraction_Ylm_mmax) cycle
      iP = iP + 1
      nodeNPots = nodeNPots + 1
      nodePairK2(iP) = k2
      nodePairL2(iP) = l2
      nodePairM(iP) = mOrb(k2) - mOrb(l2)
      nodePairPot(iP) = nodeNPots
      if (k2 .ne. l2) then
        iP = iP + 1
        nodePairK2(iP) = l2
        nodePairL2(iP) = k2
        nodePairM(iP) = -nodePairM(iP - 1)
        nodePairPot(iP) = -nodeNPots
      end if
    end do
  end do

  ! Target group of the column (iP, k): the block of the orbitals with
  ! m_k + M_iP and the body type of k (0 if there is none)
  allocate (targetGroup(nOS, nPairs))
  do iP = 1, nPairs
    do k1 = 1, nOS
      j1 = findloc(mOrb .eq. mOrb(k1) + nodePairM(iP) .and. btOrb .eq. btOrb(k1), .true., dim=1)
      targetGroup(k1, iP) = 0
      if (j1 > 0) targetGroup(k1, iP) = Orbs_blockOfOrb(j1)
    end do
  end do
  orbOfCol = spread([(k1, k1=1, nOS)], 2, nPairs)
  pairOfCol = spread([(iP, iP=1, nPairs)], 1, nOS)
  storedQ = spread(nodePairPot > 0, 1, nOS)

  allocate (nodeGroups(maxval(Orbs_blockOfOrb)))
  do iG = 1, size(nodeGroups)
    associate (g => nodeGroups(iG))
      g % orbs = pack([(k1, k1=1, nOS)], Orbs_blockOfOrb .eq. iG)
      g % colPair = [pack(pairOfCol, targetGroup .eq. iG .and. storedQ), &
                     pack(pairOfCol, targetGroup .eq. iG .and. .not. storedQ)]
      g % colOrb = [pack(orbOfCol, targetGroup .eq. iG .and. storedQ), &
                    pack(orbOfCol, targetGroup .eq. iG .and. .not. storedQ)]
      g % nColStored = count(targetGroup .eq. iG .and. storedQ)
      allocate (g % weights(size(g % colPair), size(g % orbs)))
      g % colOffset = nodeNH2Cols
      nodeNH2Cols = nodeNH2Cols + g % nColStored
      nodeMaxCols = max(nodeMaxCols, size(g % colPair))
    end associate
  end do

  ! Quadrature weights of the full measure on the node points
  nodeWeights = reshape(spread(Grid_Ylm_radialWeights(:), 2, Grid_Ylm_nNodes) * &
                        spread(Grid_Ylm_nodeWeights(:), 1, Grid_Ylm_nRadial), [Grid_Ylm_nRadial * Grid_Ylm_nNodes])

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine FillTargetGroupWeights(rdm2spatial)
  !---------------------------------------------------------------------------
  ! Weights of the correlation sum of every group of nodeGroups:
  !   weights((iP,k), j) = rdm2spatial(k, k2, j, l2),  iP = (k2, l2).
  !---------------------------------------------------------------------------
  complex(R64), intent(in), contiguous :: rdm2spatial(:, :, :, :)

  integer(I32) :: iG, iC, iJ, iP

  do iG = 1, size(nodeGroups)
    associate (g => nodeGroups(iG))
      do iJ = 1, size(g % orbs)
        do iC = 1, size(g % colPair)
          iP = g % colPair(iC)
          g % weights(iC, iJ) = rdm2spatial(g % colOrb(iC), nodePairK2(iP), g % orbs(iJ), nodePairL2(iP))
        end do
      end do
    end associate
  end do

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine FillOrbNodes(orbNodes, orbs)
  !---------------------------------------------------------------------------
  ! Orbitals on the angular nodes (channels of their declared m), flattened
  ! to (node point, orbital).
  !---------------------------------------------------------------------------
  use M_Grid_Ylm
  use M_Orbs

  complex(R64), intent(out), allocatable, target :: orbNodes(:, :)
  complex(R64), intent(in), contiguous           :: orbs(:, :)

  complex(R64), contiguous, pointer :: nodes(:, :)
  integer(I32) :: nRad, nNodes, i1

  nRad = Grid_Ylm_nRadial
  nNodes = Grid_Ylm_nNodes
  allocate (orbNodes(nRad * nNodes, Orbs_nOrbsInState))

  !$omp parallel do default(shared) private(i1, nodes)
  do i1 = 1, Orbs_nOrbsInState
    nodes(1:nRad, 1:nNodes) => orbNodes(:, i1)
    call Grid_Ylm_ToNodes(nodes, orbs(:, i1), Method_Mb_OrbBased_mOfOrb(i1))
  end do
  !$omp end parallel do

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine FillPairPotentialsYlmNodes(potNodes, orbNodes, time)
  !---------------------------------------------------------------------------
  ! Pair potentials V_(l2,k2) of the sources conjg(phi_l2) phi_k2 of the
  ! unordered pairs on the nodes, in parallel over the pairs, on the channel
  ! matrix of the single m = m_k2 - m_l2 of the pair.
  !---------------------------------------------------------------------------
  use M_Grid_Ylm
  use M_SysInteraction_Ylm

  complex(R64), intent(out), allocatable, target :: potNodes(:, :)
  complex(R64), intent(in), contiguous, target   :: orbNodes(:, :)
  real(R64), intent(in)                          :: time

  complex(R64), allocatable, target :: srcNodes(:)
  complex(R64), allocatable :: srcM(:, :), potM(:, :)
  complex(R64), contiguous, pointer :: nodes(:, :)
  integer(I32) :: nRad, nNodes, nL, iP, iL, mP
  real(R64) :: lambda

  nRad = Grid_Ylm_nRadial
  nNodes = Grid_Ylm_nNodes
  allocate (potNodes(nRad * nNodes, nodeNPots))
  ! switching factor of the interaction (S_Method_Mb_OrbBased_Switching)
  lambda = Method_Mb_OrbBased_SwitchFactor(time)

  !$omp parallel default(shared) private(iP, iL, mP, nL, srcNodes, srcM, potM, nodes)
  allocate (srcNodes(nRad * nNodes))
  !$omp do schedule(dynamic)
  do iP = 1, size(nodePairPot)
    if (nodePairPot(iP) < 0) cycle
    mP = nodePairM(iP)
    nL = SysInteraction_Ylm_lmax - abs(mP) + 1

    srcNodes(:) = conjg(orbNodes(:, nodePairL2(iP))) * orbNodes(:, nodePairK2(iP))
    nodes(1:nRad, 1:nNodes) => srcNodes
    allocate (srcM(nRad, nL), potM(nRad, nL))
    call Grid_Ylm_NodesToChannels(srcM, nodes, mP)
    do iL = 1, nL
      srcM(:, iL) = srcM(:, iL) * Grid_Ylm_radialMetricWeights(:)
    end do
    if (lambda < 1.0_R64) srcM(:, :) = lambda * srcM(:, :)

    call SysInteraction_Ylm_FillInteractionPotentialM(potM, srcM, mP, time)

    nodes(1:nRad, 1:nNodes) => potNodes(:, nodePairPot(iP))
    call Grid_Ylm_ChannelsToNodes(nodes, potM, mP)
    deallocate (srcM, potM)
  end do
  !$omp end do
  deallocate (srcNodes)
  !$omp end parallel

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine AccumulateTargetSumsYlmNodes(orbNodes, potNodes, h2Cols_, dOrbNodes_)
  !---------------------------------------------------------------------------
  ! Pointwise sums of the node path over blocks of nodes, per target group:
  ! the products X = V_iP * orbs(:,k) of the columns (iP, k) of the group
  ! enter two ZGEMMs (i, j in the group),
  !   dOrbNodes_(:,j) += X weights,
  !   h2Cols_(i, colOffset + (iP,k)) = sum w conjg(orbs(:,i)) X(:,(iP,k)),
  ! the latter for the columns of stored pairs only (FillH2FromTargetColumns).
  !---------------------------------------------------------------------------
  use M_Utils_BlasLib, only: ZGEMM
  use M_Orbs

  complex(R64), intent(in), contiguous              :: orbNodes(:, :)
  complex(R64), intent(in), contiguous              :: potNodes(:, :)
  complex(R64), intent(out), allocatable, optional  :: h2Cols_(:, :)
  complex(R64), intent(inout), contiguous, optional :: dOrbNodes_(:, :)

  !> Nodes per block of the pointwise sums (block work arrays stay in cache).
  integer(I32), parameter :: blockSize = 128

  complex(R64), allocatable :: x(:, :), orbW(:, :), dOrbBlock(:, :), h2Private(:, :)
  integer(I32) :: nOS, nPts, iB, nB, iG, iC, iJ, nC, nCS, nJ, iP, k1, j1
  logical :: h2Q, dOrbQ

  nOS = Orbs_nOrbsInState
  nPts = size(orbNodes, 1)
  h2Q = present(h2Cols_)
  dOrbQ = present(dOrbNodes_)

  if (h2Q) then
    allocate (h2Cols_(nOS, nodeNH2Cols))
    h2Cols_(:, :) = 0.0_R64
  end if

  !$omp parallel default(shared) private(iB, nB, iG, iC, iJ, nC, nCS, nJ, iP, k1, j1, x, orbW, dOrbBlock, h2Private)
  allocate (x(blockSize, nodeMaxCols))
  allocate (orbW(blockSize, nOS), dOrbBlock(blockSize, nOS))
  if (h2Q) then
    allocate (h2Private(nOS, nodeNH2Cols))
    h2Private(:, :) = 0.0_R64
  end if

  !$omp do schedule(static)
  do iB = 1, nPts, blockSize
    nB = min(blockSize, nPts - iB + 1)

    do iG = 1, size(nodeGroups)
      nC = size(nodeGroups(iG) % colPair)
      nCS = nodeGroups(iG) % nColStored
      nJ = size(nodeGroups(iG) % orbs)
      if (.not. dOrbQ) nC = nCS
      if (nC .eq. 0) cycle

      do iC = 1, nC
        iP = nodeGroups(iG) % colPair(iC)
        k1 = nodeGroups(iG) % colOrb(iC)
        if (nodePairPot(iP) > 0) then
          x(1:nB, iC) = potNodes(iB:iB + nB - 1, nodePairPot(iP)) * orbNodes(iB:iB + nB - 1, k1)
        else
          x(1:nB, iC) = conjg(potNodes(iB:iB + nB - 1, -nodePairPot(iP))) * orbNodes(iB:iB + nB - 1, k1)
        end if
      end do

      if (dOrbQ) then
        call ZGEMM('N', 'N', nB, nJ, nC, (1.0_R64, 0.0_R64), x, blockSize, nodeGroups(iG) % weights, nC, &
                   (0.0_R64, 0.0_R64), dOrbBlock, blockSize)
        do iJ = 1, nJ
          j1 = nodeGroups(iG) % orbs(iJ)
          dOrbNodes_(iB:iB + nB - 1, j1) = dOrbNodes_(iB:iB + nB - 1, j1) + dOrbBlock(1:nB, iJ)
        end do
      end if

      if (h2Q .and. nCS > 0) then
        do iJ = 1, nJ
          orbW(1:nB, iJ) = orbNodes(iB:iB + nB - 1, nodeGroups(iG) % orbs(iJ)) * nodeWeights(iB:iB + nB - 1)
        end do
        call ZGEMM('C', 'N', nJ, nCS, nB, (1.0_R64, 0.0_R64), orbW, blockSize, x, blockSize, &
                   (1.0_R64, 0.0_R64), h2Private(1, nodeGroups(iG) % colOffset + 1), nOS)
      end if
    end do
  end do
  !$omp end do

  if (h2Q) then
    !$omp critical
    h2Cols_(:, :) = h2Cols_(:, :) + h2Private(:, :)
    !$omp end critical
  end if
  !$omp end parallel

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine FillH2FromTargetColumns(h2, h2Cols)
  !---------------------------------------------------------------------------
  !   h2(i1,l2,k,k2) = 1/2 <orb_i1| V_(l2,k2) |orb_k>  for the column (iP,k)
  ! of a stored pair, the swapped pair by the Hermitian symmetry
  !   h2(k,k2,i1,l2) = conjg(h2(i1,l2,k,k2))
  ! (row k lies in the group of the swapped column (iP',i1); a diagonal pair
  ! k2 = l2 is its own swapped pair), zero elsewhere.
  ! Valid for a Hermitian metric only (precondition of the node path); a
  ! c-product (ECS) variant must compute the swapped columns instead.
  !---------------------------------------------------------------------------
  complex(R64), intent(out), allocatable :: h2(:, :, :, :)
  complex(R64), intent(in), contiguous   :: h2Cols(:, :)

  integer(I32) :: nO, iG, iC, iJ, iP, i1, k1
  complex(R64) :: h2Val

  nO = Method_Mb_OrbBased_nOrbsSum

  allocate (h2(nO, nO, nO, nO))
  h2(:, :, :, :) = 0.0_R64

  do iG = 1, size(nodeGroups)
    associate (g => nodeGroups(iG))
      do iC = 1, g % nColStored
        iP = g % colPair(iC)
        k1 = g % colOrb(iC)
        do iJ = 1, size(g % orbs)
          i1 = g % orbs(iJ)
          h2Val = 0.5_R64 * h2Cols(iJ, g % colOffset + iC)
          h2(i1, nodePairL2(iP), k1, nodePairK2(iP)) = h2Val
          if (nodePairK2(iP) .ne. nodePairL2(iP)) h2(k1, nodePairK2(iP), i1, nodePairL2(iP)) = conjg(h2Val)
        end do
      end do
    end associate
  end do

  call DuplicateRestrictedSpinBlocks(h2)

end subroutine

end submodule
