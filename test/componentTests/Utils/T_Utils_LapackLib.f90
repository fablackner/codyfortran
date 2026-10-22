!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> LapackLib_DiagonalizeHermitianBlocks on a Hermitian matrix that is block
!> diagonal in interleaved labels, with eigenvalues degenerate across blocks
!> (the situation of orbitals of different m with equal occupations or
!> norms): the eigenvectors must be exactly zero outside their block, the
!> decomposition exact to round-off and the eigenvalues ascending; with a
!> single block the result must equal LapackLib_DiagonalizeGeneric.
program T_Utils_LapackLib
  use M_Utils_Types
  use M_Utils_LapackLib
  use testdrive, only: check, error_type
  implicit none

  type(error_type), allocatable :: error

  call test_diagonalize_hermitian_blocks(error)
  if (allocated(error)) error stop 'LapackLib_DiagonalizeHermitianBlocks test failed'

contains

  subroutine test_diagonalize_hermitian_blocks(error)
    type(error_type), allocatable, intent(out) :: error

    integer(I32), parameter :: n = 6
    !> Interleaved labels: block 1 = {1, 3, 6}, block 2 = {2, 5}, block 3 = {4}
    integer(I32), parameter :: blockOf(n) = [1, 2, 1, 3, 2, 1]
    !> Eigenvalues per block: {1, 2, 3}, {1, 2}, {2} (degenerate across blocks)
    real(R64), parameter :: evalsRef(n) = [1.0_R64, 1.0_R64, 2.0_R64, 2.0_R64, 2.0_R64, 3.0_R64]

    complex(R64) :: a(n, n)
    real(R64), allocatable :: evals(:), evalsGeneric(:)
    complex(R64), allocatable :: evecs(:, :), evecsGeneric(:, :)
    integer(I32), allocatable :: idx(:)
    integer(I32) :: i, k, nOutside, iBlock
    real(R64) :: devDecomposition, devOrthonormality, devEvals
    logical :: sameGenericQ

    print *, "Testing LapackLib_DiagonalizeHermitianBlocks on a 6x6 matrix with 3 interleaved blocks..."

    a(:, :) = (0.0_R64, 0.0_R64)
    call SetBlock(a, pack([(i, i=1, n)], blockOf .eq. 1), [1.0_R64, 2.0_R64, 3.0_R64])
    call SetBlock(a, pack([(i, i=1, n)], blockOf .eq. 2), [1.0_R64, 2.0_R64])
    call SetBlock(a, pack([(i, i=1, n)], blockOf .eq. 3), [2.0_R64])

    call LapackLib_DiagonalizeHermitianBlocks(evals, evecs, a, blockOf)

    ! every eigenvector lives in exactly one block: exact zeros elsewhere
    nOutside = 0
    do k = 1, n
      iBlock = blockOf(maxloc(abs(evecs(:, k)), 1))
      do i = 1, n
        if (blockOf(i) .ne. iBlock .and. evecs(i, k) .ne. (0.0_R64, 0.0_R64)) nOutside = nOutside + 1
      end do
    end do

    devEvals = maxval(abs(evals - evalsRef))
    devDecomposition = maxval(abs(matmul(a, evecs) - evecs * spread(evals, 1, n)))
    devOrthonormality = maxval(abs(matmul(conjg(transpose(evecs)), evecs) - Identity(n)))

    print *, "nonzero entries outside the block: ", nOutside
    print *, "eigenvalue deviation:              ", devEvals
    print *, "decomposition residual:            ", devDecomposition
    print *, "orthonormality deviation:          ", devOrthonormality

    call check(error, nOutside, 0)
    if (allocated(error)) return
    call check(error, devEvals < 1.0e-13_R64)
    if (allocated(error)) return
    call check(error, devDecomposition < 1.0e-13_R64)
    if (allocated(error)) return
    call check(error, devOrthonormality < 1.0e-13_R64)
    if (allocated(error)) return

    ! a single block is the plain diagonalization (bitwise)
    idx = pack([(i, i=1, n)], blockOf .eq. 1)
    call LapackLib_DiagonalizeHermitianBlocks(evals, evecs, a(idx, idx), [(1, i=1, size(idx))])
    call LapackLib_DiagonalizeGeneric(evalsGeneric, evecsGeneric, a(idx, idx), .true.)
    sameGenericQ = all(evals .eq. evalsGeneric) .and. all(evecs .eq. evecsGeneric)

    print *, "single block equals DiagonalizeGeneric: ", sameGenericQ

    call check(error, sameGenericQ)

  end subroutine test_diagonalize_hermitian_blocks

  !> a(idx, idx) = H diag(evals) H with the Householder reflection
  !> H = 1 - 2 v v^H / (v^H v) (unitary and Hermitian, mixes all entries).
  subroutine SetBlock(a, idx, evals)
    complex(R64), intent(inout) :: a(:, :)
    integer(I32), intent(in) :: idx(:)
    real(R64), intent(in) :: evals(:)

    complex(R64), allocatable :: v(:), h(:, :)
    integer(I32) :: nB, i

    nB = size(idx)
    allocate (v(nB))
    do i = 1, nB
      v(i) = cmplx(1.0_R64 + 0.3_R64 * i, 0.7_R64 - 0.2_R64 * i**2, kind=R64)
    end do
    h = Identity(nB) - 2.0_R64 * spread(v, 2, nB) * spread(conjg(v), 1, nB) / real(dot_product(v, v), kind=R64)

    a(idx, idx) = matmul(h, matmul(Diagonal(evals), h))

  end subroutine

  function Identity(n) result(res)
    integer(I32), intent(in) :: n
    complex(R64) :: res(n, n)

    integer(I32) :: i

    res(:, :) = (0.0_R64, 0.0_R64)
    do i = 1, n
      res(i, i) = (1.0_R64, 0.0_R64)
    end do

  end function

  function Diagonal(d) result(res)
    real(R64), intent(in) :: d(:)
    complex(R64) :: res(size(d), size(d))

    integer(I32) :: i

    res(:, :) = (0.0_R64, 0.0_R64)
    do i = 1, size(d)
      res(i, i) = cmplx(d(i), 0.0_R64, kind=R64)
    end do

  end function

end program T_Utils_LapackLib
