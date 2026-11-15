! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The closure of the scheme plus the connected triples; see
!> A_Reconstruction_Singlet_Triples.md.
!>
!> Steps of Reconstruction_Singlet_Triples_FillCumul3:
!>   1. semicanonical orbitals: the nBodies(1) most occupied natural orbitals
!>      of rdm1a and the others, each set rotated to diagonalize its block of
!>      the Fock operator f = h1 + J - K of rdm1a;
!>   2. 2-cumulant and interaction in these orbitals;
!>   3. source of the up-up-down block (FillSource);
!>   4. amplitudes of the variant;
!>   5. the block and its adjoint, back in the orbitals of the 2-RDM.
submodule(M_Reconstruction_Singlet_Triples) S_Reconstruction_Singlet_Triples
  use M_Reconstruction_Singlet, only: I_Reconstruction_Singlet_FillRdm3aab

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> The closure that this extension wraps.
  procedure(I_Reconstruction_Singlet_FillRdm3aab), pointer :: InnerFillRdm3aab => null()

  !> The Hamiltonian of Reconstruction_SetHamiltonian in spatial orbitals:
  !> h1 of the up block and the interaction <pq|rs> (h2 = <pq|rs>/2).
  complex(R64), allocatable :: spatialH1(:, :), spatialW(:, :, :, :)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Singlet_Triples_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Reconstruction
    use M_Reconstruction_Singlet
    use M_Reconstruction_Singlet_Triples_MoellerPlesset

    call Say_Fabricate("reconstruction.singlet.triples")

    Reconstruction_Singlet_Triples_FillAmplitudes => null()
    if (Json_GetExistence("reconstruction.singlet.triples.moellerPlesset")) then
      call Reconstruction_Singlet_Triples_MoellerPlesset_Fabricate
    else
      error stop "reconstruction.singlet.triples is missing one of: moellerPlesset"
    end if
    if (.not. associated(Reconstruction_Singlet_Triples_FillAmplitudes)) &
      error stop "reconstruction.singlet.triples: the variant binds no amplitudes"

    if (allocated(spatialH1)) deallocate (spatialH1)
    if (allocated(spatialW)) deallocate (spatialW)

    InnerFillRdm3aab => Reconstruction_Singlet_FillRdm3aab
    Reconstruction_Singlet_FillRdm3aab => FillRdm3aab
    Reconstruction_SetHamiltonian => SetHamiltonian

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Keeps the up block of h1 and the spatial interaction (the up-down block
  !> of h2, the one the collision operator contracts).
  subroutine SetHamiltonian(h1, h2)
    complex(R64), intent(in), contiguous :: h1(:, :)
    complex(R64), intent(in), contiguous :: h2(:, :, :, :)

    integer(I32) :: n

    n = size(h1, 1) / 2
    spatialH1 = h1(1:n, 1:n)
    spatialW = 2.0_R64 * h2(1:n, n + 1:, 1:n, n + 1:)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm3aab(rdm3aab, rdm1a, rdm2ab)
    complex(R64), intent(out), allocatable :: rdm3aab(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm1a(:, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    complex(R64), allocatable :: cumul3(:, :, :, :, :, :)

    call InnerFillRdm3aab(rdm3aab, rdm1a, rdm2ab)
    call Reconstruction_Singlet_Triples_FillCumul3(cumul3, rdm1a, rdm2ab)
    rdm3aab(:, :, :, :, :, :) = rdm3aab(:, :, :, :, :, :) + cumul3(:, :, :, :, :, :)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Singlet_Triples_FillCumul3(cumul3, rdm1a, rdm2ab)
    use M_Utils_RdmBasisChange
    use M_Utils_RdmCumulants
    use M_Method_Mb

    complex(R64), intent(out), allocatable :: cumul3(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm1a(:, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    complex(R64), allocatable :: orbs(:, :), toOrbs(:, :), rdm1(:, :), rdm2(:, :, :, :), cumul2(:, :, :, :), w(:, :, :, :)
    complex(R64), allocatable :: source(:, :, :, :, :, :), amplitudes(:, :, :, :, :, :), block(:, :, :, :, :, :)
    real(R64), allocatable :: epsOcc(:), epsVirt(:)
    integer(I32) :: n, nOcc, nVirt, i1, i2, i3, j1, j2, j3

    n = size(rdm1a, 1)
    nOcc = Method_Mb_nBodies(1)
    nVirt = n - nOcc
    if (nOcc < 1 .or. nVirt < 1) then
      allocate (cumul3(n, n, n, n, n, n), source=cmplx(0.0_R64, 0.0_R64, R64))
      return
    end if
    if (.not. allocated(spatialH1)) &
      error stop "reconstruction.singlet.triples: no Hamiltonian; the method must call Reconstruction_SetHamiltonian "// &
      "before the closure (td2rdmStaticOrbs does)"
    if (size(spatialH1, 1) .ne. n) error stop "reconstruction.singlet.triples: the Hamiltonian has another number of orbitals"

    call FillSemicanonicalOrbitals(orbs, epsOcc, epsVirt, rdm1a, nOcc)

    ! <pq|rs> transforms as a ket-first RDM does: conjugated orbitals on p, q
    toOrbs = conjg(transpose(orbs))
    rdm1 = matmul(toOrbs, matmul(rdm1a, orbs))
    rdm2 = rdm2ab
    call RdmBasisChange_Transform(rdm2, 4, toOrbs)
    call RdmCumulants_FillCumul2ab(cumul2, rdm1, rdm1, rdm2)
    w = spatialW
    call RdmBasisChange_Transform(w, 4, toOrbs)

    call FillSource(source, cumul2, w, nOcc)
    allocate (amplitudes(nOcc, nOcc, nOcc, nVirt, nVirt, nVirt))
    call Reconstruction_Singlet_Triples_FillAmplitudes(amplitudes, source, epsOcc, epsVirt)

    ! the block ooo|vvv in the orbitals of the 2-RDM plus its adjoint
    allocate (block(n, n, n, n, n, n), cumul3(n, n, n, n, n, n))
    call RdmBasisChange_EmbedBlock(block, amplitudes, 6, orbs(:, :nOcc), orbs(:, nOcc + 1:))
    !$omp parallel do collapse(3) private(i1, i2, i3, j1, j2, j3)
    do j3 = 1, n
      do j2 = 1, n
        do j1 = 1, n
          do i3 = 1, n
            do i2 = 1, n
              do i1 = 1, n
                cumul3(i1, i2, i3, j1, j2, j3) = block(i1, i2, i3, j1, j2, j3) + conjg(block(j1, j2, j3, i1, i2, i3))
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Columns of orbs: the nOcc most occupied natural orbitals of rdm1a, then
  !> the others, each set rotated to diagonalize its block of the Fock
  !> operator f(i,j) = h1(i,j) + sum_kl (2 <ik|jl> - <ik|lj>) rdm1a(l,k);
  !> epsOcc and epsVirt its diagonal.
  subroutine FillSemicanonicalOrbitals(orbs, epsOcc, epsVirt, rdm1a, nOcc)
    use M_Utils_LapackLib
    use M_Utils_RdmDiagonalize

    complex(R64), intent(out), allocatable :: orbs(:, :)
    real(R64), intent(out), allocatable    :: epsOcc(:)
    real(R64), intent(out), allocatable    :: epsVirt(:)
    complex(R64), intent(in), contiguous   :: rdm1a(:, :)
    integer(I32), intent(in)               :: nOcc

    complex(R64), allocatable :: fock(:, :), natOrbs(:, :), evecs(:, :)
    real(R64), allocatable :: natOcc(:)
    integer(I32) :: n, i, j, k, l

    n = size(rdm1a, 1)
    fock = spatialH1
    !$omp parallel do collapse(2) private(i, j, k, l)
    do j = 1, n
      do i = 1, n
        do l = 1, n
          do k = 1, n
            fock(i, j) = fock(i, j) + (2.0_R64 * spatialW(i, k, j, l) - spatialW(i, k, l, j)) * rdm1a(l, k)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call RdmDiagonalize_Rdm1(natOcc, natOrbs, rdm1a)
    fock = matmul(conjg(transpose(natOrbs)), matmul(fock, natOrbs))
    allocate (orbs(n, n))
    call LapackLib_DiagonalizeGeneric(epsOcc, evecs, fock(:nOcc, :nOcc), .true.)
    orbs(:, :nOcc) = matmul(natOrbs(:, :nOcc), evecs)
    call LapackLib_DiagonalizeGeneric(epsVirt, evecs, fock(nOcc + 1:, nOcc + 1:), .true.)
    orbs(:, nOcc + 1:) = matmul(natOrbs(:, nOcc + 1:), evecs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> The up-up-down block of the spin-orbital source
  !>   P(k/ij) P(a/bc) sum_e L(ij,ae) <ek||bc> - P(i/jk) P(c/ab) sum_m L(im,ab) <jk||mc>,
  !> P(x/yz) f = f - f(x <-> y) - f(x <-> z): the second-order connected
  !> triples with the first-order doubles replaced by the 2-cumulant L. With
  !> the up-down cumulant cumul2 of a singlet (cumul2(p,q,r,s) =
  !> cumul2(q,p,s,r)) and w(p,q,r,s) = <pq|rs>, its 18 terms are entries of
  !> two matrix products,
  !>   pLine((p,q,r),(x,y,z)) = sum_e cumul2(p,q,r,e) w(e,x,y,z)   particle line,
  !>   hLine((p,r,s),(x,y,z)) = sum_m cumul2(p,m,r,s) w(x,y,m,z)   hole line
  !> (e unoccupied, m occupied): 12 entries of each.
  subroutine FillSource(source, cumul2, w, nOcc)
    use M_Utils_BlasLib
    use M_Utils_Constants, only: CONE, CZERO

    complex(R64), intent(out), allocatable :: source(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous   :: cumul2(:, :, :, :)
    complex(R64), intent(in), contiguous   :: w(:, :, :, :)
    integer(I32), intent(in)               :: nOcc

    complex(R64), allocatable :: pLine(:, :, :, :, :, :), hLine(:, :, :, :, :, :)
    complex(R64), allocatable :: doubles(:, :, :, :), wLine(:, :, :, :)
    integer(I32) :: nVirt, i, j, k, a, b, c

    nVirt = size(w, 1) - nOcc
    allocate (pLine(nOcc, nOcc, nVirt, nOcc, nVirt, nVirt), hLine(nOcc, nVirt, nVirt, nOcc, nOcc, nVirt))
    allocate (source(nOcc, nOcc, nOcc, nVirt, nVirt, nVirt))

    ! doubles(p,q,r,e) wLine(e,x,y,z)
    doubles = cumul2(:nOcc, :nOcc, nOcc + 1:, nOcc + 1:)
    wLine = w(nOcc + 1:, :nOcc, nOcc + 1:, nOcc + 1:)
    call ZGEMM('N', 'N', nOcc**2 * nVirt, nOcc * nVirt**2, nVirt, CONE, doubles, nOcc**2 * nVirt, wLine, nVirt, &
               CZERO, pLine, nOcc**2 * nVirt)
    ! doubles(p,r,s,m) wLine(m,x,y,z)
    doubles = reshape(cumul2(:nOcc, :nOcc, nOcc + 1:, nOcc + 1:), [nOcc, nVirt, nVirt, nOcc], order=[1, 4, 2, 3])
    wLine = reshape(w(:nOcc, :nOcc, :nOcc, nOcc + 1:), [nOcc, nOcc, nOcc, nVirt], order=[2, 3, 1, 4])
    call ZGEMM('N', 'N', nOcc * nVirt**2, nOcc**2 * nVirt, nOcc, CONE, doubles, nOcc * nVirt**2, wLine, nOcc, &
               CZERO, hLine, nOcc * nVirt**2)

    !$omp parallel do collapse(2) private(i, j, k, a, b, c)
    do c = 1, nVirt
      do b = 1, nVirt
        do a = 1, nVirt
          do k = 1, nOcc
            do j = 1, nOcc
              do i = 1, nOcc
                source(i, j, k, a, b, c) = &
                  pLine(i, j, a, k, b, c) - pLine(j, i, a, k, b, c) - pLine(i, j, b, k, a, c) + pLine(j, i, b, k, a, c) &
                  - pLine(j, k, a, i, c, b) + pLine(j, k, b, i, c, a) + pLine(i, k, a, j, c, b) - pLine(i, k, b, j, c, a) &
                  + pLine(k, j, c, i, b, a) - pLine(k, j, c, i, a, b) - pLine(k, i, c, j, b, a) + pLine(k, i, c, j, a, b) &
                  - hLine(i, a, b, j, k, c) + hLine(i, b, a, j, k, c) + hLine(j, a, b, i, k, c) - hLine(j, b, a, i, k, c) &
                  + hLine(i, b, c, k, j, a) - hLine(i, a, c, k, j, b) - hLine(j, b, c, k, i, a) + hLine(j, a, c, k, i, b) &
                  - hLine(k, c, b, j, i, a) + hLine(k, c, b, i, j, a) + hLine(k, c, a, j, i, b) - hLine(k, c, a, i, j, b)
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine

end submodule
