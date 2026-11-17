!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Random Hamiltonians and dense CI matrices, shared by the component tests
! of coeffs.hubbard (T_Coeffs_Hubbard_VsGeneric, T_Coeffs_Hubbard_SpinSectors);
! included after `contains`.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Random Hermitian n x n matrix, entries uniform in the unit square around 0.
subroutine FillRandomHermitian(a, n)
  use M_Utils_Types

  complex(R64), intent(out), allocatable :: a(:, :)
  integer(I32), intent(in) :: n

  real(R64) :: re(n, n), im(n, n)

  call random_number(re)
  call random_number(im)
  a = cmplx(re - 0.5_R64, im - 0.5_R64, kind=R64)
  a = 0.5_R64 * (a + conjg(transpose(a)))

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Spin-orbital h1 with the blocks hUp and hDn of the two spins and the
!> spin-orbital h2 of the on-site interaction sum_s u(s) n_{s,up} n_{s,dn}:
!> h2(s up, s dn, s up, s dn) = h2(s dn, s up, s dn, s up) = u(s) / 2.
subroutine FillLatticeHamiltonian(h1, h2, hUp, hDn, u)
  use M_Utils_Types
  use M_Method_Mb_OrbBased

  complex(R64), intent(out), allocatable :: h1(:, :)
  complex(R64), intent(out), allocatable :: h2(:, :, :, :)
  complex(R64), intent(in), contiguous :: hUp(:, :)
  complex(R64), intent(in), contiguous :: hDn(:, :)
  real(R64), intent(in), contiguous :: u(:)

  integer(I32) :: nO, n, oUp, oDn, s

  nO = Method_Mb_OrbBased_nOrbsSum
  n = size(u)
  oUp = Method_Mb_OrbBased_nOrbsStart(1) - 1
  oDn = Method_Mb_OrbBased_nOrbsStart(2) - 1

  allocate (h1(nO, nO), h2(nO, nO, nO, nO))
  h1(:, :) = 0.0_R64
  h1(oUp + 1:oUp + n, oUp + 1:oUp + n) = hUp
  h1(oDn + 1:oDn + n, oDn + 1:oDn + n) = hDn

  h2(:, :, :, :) = 0.0_R64
  do s = 1, n
    h2(oUp + s, oDn + s, oUp + s, oDn + s) = 0.5_R64 * u(s)
    h2(oDn + s, oUp + s, oDn + s, oUp + s) = 0.5_R64 * u(s)
  end do

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Dense CI matrix of H1 + H2 of the bound Coeffs backend, column j the
!> Hamiltonian applied to the unit vector j.
subroutine FillDenseHamiltonian(hMat, h1, h2)
  use M_Utils_Types
  use M_Coeffs

  complex(R64), intent(out), allocatable :: hMat(:, :)
  complex(R64), intent(in), contiguous :: h1(:, :)
  complex(R64), intent(in), contiguous :: h2(:, :, :, :)

  complex(R64), allocatable :: unit(:), column(:)
  integer(I32) :: j

  allocate (hMat(Coeffs_nCoeffs, Coeffs_nCoeffs), unit(Coeffs_nCoeffs), column(Coeffs_nCoeffs))
  do j = 1, Coeffs_nCoeffs
    unit(:) = 0.0_R64
    unit(j) = 1.0_R64
    column(:) = 0.0_R64
    call Coeffs_ApplyH1FillRdm1(unit, dCoeffs_=column, h1_=h1)
    call Coeffs_ApplyH2FillRdm2(unit, dCoeffs_=column, h2_=h2)
    hMat(:, j) = column
  end do

end subroutine
