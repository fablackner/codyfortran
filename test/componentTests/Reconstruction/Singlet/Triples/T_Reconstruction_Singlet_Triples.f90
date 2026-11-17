!> Connected triples of the singlet closure (reconstruction.singlet.triples)
!> on a model with 5 orbitals and 2 electrons per spin (complex Hermitian h1,
!> interaction and 2-cumulant from closed formulas):
!>   - the added 3-cumulant T (Reconstruction_Singlet_Triples_FillCumul3) is
!>     antisymmetric, Hermitian and has no onefold contraction (so ccQ leaves
!>     it unchanged);
!>   - covariance: with the 2-RDM and the Hamiltonian in rotated orbitals
!>     (complex unitary) T is the rotated T, which fails for a wrong
!>     conjugation or index order in the source;
!>   - regression: |T|, pinned from the spin-orbital formula of the source.
program T_Reconstruction_Singlet_Triples
  use M_Utils_Types
  use M_Utils_Json
  use M_Utils_LapackLib
  use M_Utils_RdmBasisChange
  use M_Reconstruction
  use M_Reconstruction_Singlet_Triples
  use M_Method_Mb
  use testdrive, only: check, error_type
  implicit none

  integer(I32), parameter :: n = 5, nUp = 2
  character(len=*), parameter :: dir = "test/componentTests/Reconstruction/Singlet/Triples/"
  character(len=*), parameter :: variants(1) = [character(len=14) :: "MoellerPlesset"]
  real(R64), parameter :: pinnedNorms(1) = [2.7035148855974672e-1_R64]

  type(error_type), allocatable :: error
  complex(R64), allocatable :: rdm1a(:, :), rdm2ab(:, :, :, :), h1(:, :), w(:, :, :, :), rot(:, :)
  complex(R64), allocatable :: t(:, :, :, :, :, :), tRot(:, :, :, :, :, :)
  real(R64) :: tNorm
  integer(I32) :: iVariant

  Method_Mb_nBodyTypes = 2
  allocate (Method_Mb_nBodies(2))
  Method_Mb_nBodies = [nUp, nUp]

  call FillModel(rdm1a, rdm2ab, h1, w, rot)

  do iVariant = 1, size(variants)
    call Json_LoadJsonFile(manualFileName_=dir//"T_Reconstruction_Singlet_Triples_"//trim(variants(iVariant))//".json", &
                           relativeToProjectDirQ_=.true.)
    call Reconstruction_Fabricate

    call FillTriples(t, rdm1a, rdm2ab, h1, w)
    tNorm = norm2(abs(t))
    write (*, '(A,A,ES24.16)') trim(variants(iVariant)), ": |T| = ", tNorm

    call check(error, tNorm > 1e-4_R64, "the triples vanish")
    if (.not. allocated(error)) call check(error, maxval(abs(t + reshape(t, shape(t), order=[2, 1, 3, 4, 5, 6]))) < 1e-13_R64, &
                                           "not antisymmetric in (i1, i2)")
    if (.not. allocated(error)) call check(error, maxval(abs(t + reshape(t, shape(t), order=[1, 2, 3, 5, 4, 6]))) < 1e-13_R64, &
                                           "not antisymmetric in (j1, j2)")
    if (.not. allocated(error)) &
      call check(error, maxval(abs(t - conjg(reshape(t, shape(t), order=[4, 5, 6, 1, 2, 3])))) < 1e-13_R64, "not Hermitian")
    if (.not. allocated(error)) call check(error, MaxContraction(t) < 1e-13_R64, "onefold contractions do not vanish")

    ! rotated orbitals: RDMs and integrals transform alike, h1 -> rot h1 rot^+
    if (.not. allocated(error)) then
      block
        complex(R64), allocatable :: rdm1aRot(:, :), rdm2abRot(:, :, :, :), h1Rot(:, :), wRot(:, :, :, :)
        rdm1aRot = matmul(rot, matmul(rdm1a, conjg(transpose(rot))))
        h1Rot = matmul(rot, matmul(h1, conjg(transpose(rot))))
        rdm2abRot = rdm2ab
        wRot = w
        call RdmBasisChange_Transform(rdm2abRot, 4, rot)
        call RdmBasisChange_Transform(wRot, 4, rot)
        call FillTriples(tRot, rdm1aRot, rdm2abRot, h1Rot, wRot)
        call RdmBasisChange_Transform(t, 6, rot)
        write (*, '(A,ES10.2)') "  covariance error ", maxval(abs(tRot - t))
        call check(error, maxval(abs(tRot - t)) < 1e-9_R64 * tNorm, "not covariant under an orbital rotation")
      end block
    end if

    if (.not. allocated(error)) call check(error, abs(tNorm - pinnedNorms(iVariant)) < 1e-9_R64 * tNorm, "|T| changed")
    if (allocated(error)) then
      write (*, '(A)') trim(variants(iVariant))//": "//error % message
      error stop "T_Reconstruction_Singlet_Triples failed"
    end if
  end do

contains

  !> The triples cumulant for the Hamiltonian h1, w = <pq|rs> (spatial
  !> orbitals).
  subroutine FillTriples(t, rdm1a, rdm2ab, h1, w)
    complex(R64), intent(out), allocatable :: t(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm1a(:, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)
    complex(R64), intent(in), contiguous   :: h1(:, :)
    complex(R64), intent(in), contiguous   :: w(:, :, :, :)

    complex(R64), allocatable :: h1So(:, :), h2So(:, :, :, :)

    ! spin orbitals, up block first; h2 = <pq|rs>/2 in the four spin blocks
    allocate (h1So(2 * n, 2 * n), h2So(2 * n, 2 * n, 2 * n, 2 * n))
    h1So = 0.0_R64
    h1So(:n, :n) = h1
    h1So(n + 1:, n + 1:) = h1
    h2So = 0.0_R64
    h2So(:n, :n, :n, :n) = 0.5_R64 * w
    h2So(n + 1:, :n, n + 1:, :n) = 0.5_R64 * w
    h2So(:n, n + 1:, :n, n + 1:) = 0.5_R64 * w
    h2So(n + 1:, n + 1:, n + 1:, n + 1:) = 0.5_R64 * w

    call Reconstruction_SetHamiltonian(h1So, h2So)
    call Reconstruction_Singlet_Triples_FillCumul3(t, rdm1a, rdm2ab)

  end subroutine

  !> Largest element of the four onefold contractions of an M3aab tensor.
  real(R64) function MaxContraction(t)
    complex(R64), intent(in) :: t(:, :, :, :, :, :)

    complex(R64), allocatable :: c(:, :, :, :, :)
    integer(I32) :: k

    allocate (c(n, n, n, n, 4))
    c = 0.0_R64
    do k = 1, n
      c(:, :, :, :, 1) = c(:, :, :, :, 1) + t(k, :, :, k, :, :)
      c(:, :, :, :, 2) = c(:, :, :, :, 2) + t(:, k, :, :, :, k)
      c(:, :, :, :, 3) = c(:, :, :, :, 3) + t(:, :, k, :, k, :)
      c(:, :, :, :, 4) = c(:, :, :, :, 4) + t(:, :, k, :, :, k)
    end do
    MaxContraction = maxval(abs(c))

  end function

  !> A deterministic model: Hermitian h1 with a gap, w(p,q,r,s) = <pq|rs>
  !> (Hermitian, symmetric under the exchange of the particles), a singlet
  !> rdm2ab = rdm1a x rdm1a + cumulant with two strongly occupied natural
  !> orbitals, and a complex unitary rot.
  subroutine FillModel(rdm1a, rdm2ab, h1, w, rot)
    complex(R64), intent(out), allocatable :: rdm1a(:, :), rdm2ab(:, :, :, :), h1(:, :), w(:, :, :, :), rot(:, :)

    complex(R64), allocatable :: a(:, :), c(:, :, :, :)
    real(R64), allocatable :: evals(:)
    real(R64), parameter :: occupations(n) = [0.97_R64, 0.94_R64, 0.05_R64, 0.03_R64, 0.01_R64]
    integer(I32) :: p, q, r, s

    allocate (a(n, n), h1(n, n), rdm1a(n, n), w(n, n, n, n), c(n, n, n, n), rdm2ab(n, n, n, n))

    do q = 1, n
      do p = 1, n
        a(p, q) = cmplx(Noise(p + q, p * q), merge(1, -1, p < q) * Noise(p * q + 1, p + q) * merge(0, 1, p .eq. q), R64)
      end do
    end do
    a = 0.5_R64 * (a + conjg(transpose(a)))
    call LapackLib_DiagonalizeGeneric(evals, rot, a, .true.)

    ! h1: levels -2, -1.5, 0.5, 1, 1.5 in orbitals close to the natural ones
    h1 = 0.05_R64 * a
    do p = 1, n
      h1(p, p) = h1(p, p) + merge(-2.5_R64, -1.0_R64, p <= nUp) + 0.5_R64 * p
      rdm1a(:, p) = 0.0_R64
      rdm1a(p, p) = occupations(p)
    end do
    rdm1a = rdm1a + 0.01_R64 * a

    do s = 1, n
      do r = 1, n
        do q = 1, n
          do p = 1, n
            w(p, q, r, s) = 0.3_R64 * cmplx(Noise(p + r, q + s) + Noise(q + s, p + r), &
                                            Noise(p + 2 * q, r + 2 * s) - Noise(r + 2 * s, p + 2 * q), R64)
            c(p, q, r, s) = 0.02_R64 * cmplx(Noise(p * r, q * s) + Noise(q * s, p * r), &
                                             Noise(p + 3 * q, r + 3 * s) - Noise(r + 3 * s, p + 3 * q), R64)
          end do
        end do
      end do
    end do
    ! Hermitian and symmetric under (p,r) <-> (q,s)
    w = 0.5_R64 * (w + reshape(w, shape(w), order=[2, 1, 4, 3]))
    w = 0.5_R64 * (w + conjg(reshape(w, shape(w), order=[3, 4, 1, 2])))
    c = 0.5_R64 * (c + reshape(c, shape(c), order=[2, 1, 4, 3]))
    c = 0.5_R64 * (c + conjg(reshape(c, shape(c), order=[3, 4, 1, 2])))

    do s = 1, n
      do r = 1, n
        do q = 1, n
          do p = 1, n
            rdm2ab(p, q, r, s) = rdm1a(p, r) * rdm1a(q, s) + c(p, q, r, s)
          end do
        end do
      end do
    end do

  end subroutine

  !> A fixed number in (-1, 1) for a pair of integers.
  pure real(R64) function Noise(i, j)
    integer(I32), intent(in) :: i, j

    Noise = sin(12.9898_R64 * i + 78.233_R64 * j + 0.37_R64 * i * j)

  end function

end program
