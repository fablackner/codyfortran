!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Round trip of Method_Mb_OrbBased_OrthonormalizeKeepingState, shared by
! T_He1d_09 (unrestricted) and T_He1d_10 (restricted); included after
! `contains`.
!
! The orthonormal orbitals phi and CI vector c of a full CI MCTDHX state are
! rewritten as phi M (M Hermitian positive, block diagonal in Orbs_blockOfOrb)
! with c' = exp(-(log M)^) c, the same many-body state. Its overlap matrix is
! M^2, so the Loewdin step must return phi and c.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine CheckOrthonormalizeKeepingState(jsonFileName)
  use M_Utils_Json
  use M_Utils_Say
  use M_Utils_LapackLib
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  character(len=*), intent(in) :: jsonFileName

  type(error_type), allocatable :: error
  complex(R64), allocatable :: orbsRef(:, :), coeffsRef(:), logM(:, :), kFull(:, :), evecs(:, :), term(:), termNew(:)
  real(R64), allocatable :: evals(:)
  integer(I32) :: nOS, i, j, k

  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup

  ! orthonormal reference state
  call Method_Mb_OrbBased_OrthonormalizeKeepingState(Method_state)
  call Coeffs_Normalize(Coeffs_coeffs)
  orbsRef = Orbs_orbs
  coeffsRef = Coeffs_coeffs
  nOS = Orbs_nOrbsInState

  ! log M: Hermitian, block diagonal, ||log M|| ~ 0.05
  allocate (logM(nOS, nOS))
  logM(:, :) = 0.0_R64
  do j = 1, nOS
    do i = 1, nOS
      if (Orbs_blockOfOrb(i) .ne. Orbs_blockOfOrb(j)) cycle
      logM(i, j) = cmplx(0.01_R64 * (i + j), 0.005_R64 * (i - j), kind=R64)
    end do
  end do
  call LapackLib_DiagonalizeHermitianBlocks(evals, evecs, logM, Orbs_blockOfOrb)

  ! phi <- phi M
  Orbs_orbs(:, :) = matmul(orbsRef, matmul(evecs, spread(exp(evals), 2, nOS) * conjg(transpose(evecs))))

  ! c <- exp(-(log M)^) c in the spin-orbital basis
  allocate (kFull(Method_Mb_OrbBased_nOrbsSum, Method_Mb_OrbBased_nOrbsSum))
  kFull(:, :) = 0.0_R64
  kFull(1:nOS, 1:nOS) = -logM
  if (Orbs_restrictedQ) kFull(nOS + 1:, nOS + 1:) = -logM
  term = Coeffs_coeffs
  allocate (termNew, mold=term)
  do k = 1, 60
    termNew(:) = 0.0_R64
    call Coeffs_ApplyH1FillRdm1(term, dCoeffs_=termNew, h1_=kFull)
    term = termNew / k
    Coeffs_coeffs = Coeffs_coeffs + term
    if (maxval(abs(term)) < 1.0e-16_R64) exit
  end do

  call Method_Mb_OrbBased_OrthonormalizeKeepingState(Method_state)

  print '(A,ES10.2)', " max |phi - phi_ref| = ", maxval(abs(Orbs_orbs - orbsRef))
  print '(A,ES10.2)', " max |c - c_ref|     = ", maxval(abs(Coeffs_coeffs - coeffsRef))

  call check(error, maxval(abs(Orbs_orbs - orbsRef)) < 1e-12_R64)
  if (allocated(error)) error stop "OrthonormalizeKeepingState: orbitals not recovered"
  call check(error, maxval(abs(Coeffs_coeffs - coeffsRef)) < 1e-12_R64)
  if (allocated(error)) error stop "OrthonormalizeKeepingState: CI vector not recovered"

end subroutine
