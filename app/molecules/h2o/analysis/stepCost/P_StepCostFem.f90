!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Cost of one time derivative of MCTDHF (method mctdhx) or TD-2RDM (method
!> td2rdm) for H2O on grid.fem, as a function of the number of orbitals:
!> the measurement behind the crossover table of app/README.md.
!>
!> The orbitals are those of orbsInit (STOs; orthonormalized here), the
!> many-body part a configuration (unary CI vector) or the corresponding
!> Slater-determinant 2-RDM (twoRdmInit.excited): the cost does not depend on
!> the values. Timed, program.nRepeat times each (wall clock, average):
!>   total      Method_TimeDerivative
!>   CI         Method_Mb_OrbBased_ApplyCiHamiltonian (mctdhx)
!>   2-RDM      TwoRdm_ApplyCommutatorH1/H2 + Reconstruction_AddCollisionOp
!>              with the h1, h2 of the last derivative (td2rdm)
!> The rest of total is the grid part (one-body operator, pair potentials by
!> Poisson solves, mean fields, projection), the same for both methods.
program P_StepCostFem
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_TwoRdm
  use M_TwoRdmInit
  use M_Reconstruction
  use M_Method
  use M_Method_Mb_OrbBased

  implicit none

  logical      :: ciQ
  integer(I32) :: nRepeat, r
  integer(I64) :: c0, c1, cRate
  real(R64)    :: tTotal, tPart
  complex(R64), allocatable :: dState(:), dCoeffs(:), dTwo(:)

  call Say_Hello
  call Json_LoadJsonFile()
  ciQ = Json_GetExistence("method.mb.orbBased.mctdhx")
  nRepeat = Json_Get("program.nRepeat", 2)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  if (ciQ) then
    call ConfigList_Fabricate
    call Coeffs_Fabricate
    call CoeffsInit_Fabricate
  else
    call TwoRdm_Fabricate
    call TwoRdmInit_Fabricate
    call Reconstruction_Fabricate
  end if

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  if (ciQ) then
    call ConfigList_Setup
    call Coeffs_Setup
    call CoeffsInit_Setup
  else
    call TwoRdm_Setup
    call TwoRdmInit_Setup
  end if
  call Method_Setup
  if (.not. ciQ) call Reconstruction_Setup
  call Orbs_Orthonormalize(Orbs_orbs)

  allocate (dState, mold=Method_state)
  call system_clock(c0, cRate)
  do r = 1, nRepeat
    call Method_TimeDerivative(dState, Method_state, 0.0_R64)
  end do
  call system_clock(c1)
  tTotal = real(c1 - c0, R64) / cRate / nRepeat

  ! h1, h2 of the orbitals for the many-body part alone
  call Method_Mb_OrbBased_FillH1(Method_Mb_OrbBased_h1, Orbs_orbs, 0.0_R64)
  call Method_Mb_OrbBased_FillH2(Method_Mb_OrbBased_h2, Orbs_orbs, 0.0_R64)

  call system_clock(c0)
  if (ciQ) then
    allocate (dCoeffs(Coeffs_nCoeffs))
    do r = 1, nRepeat
      call Method_Mb_OrbBased_ApplyCiHamiltonian(dCoeffs, Coeffs_coeffs, 0.0_R64)
    end do
  else
    allocate (dTwo(TwoRdm_nElements))
    do r = 1, nRepeat
      dTwo = 0.0_R64
      call TwoRdm_ApplyCommutatorH1(dTwo, TwoRdm_twoRdmState, Method_Mb_OrbBased_h1)
      call TwoRdm_ApplyCommutatorH2(dTwo, TwoRdm_twoRdmState, Method_Mb_OrbBased_h2)
      call Reconstruction_AddCollisionOp(dTwo, TwoRdm_twoRdmState, Method_Mb_OrbBased_h2)
    end do
  end if
  call system_clock(c1)
  tPart = real(c1 - c0, R64) / cRate / nRepeat

  write (*, '(A,A,A,I0,A,I0,A,I0,A,F10.3,A,F10.3,A,F10.3)') " STEPCOST method=", merge("mctdhf", "td2rdm", ciQ), &
    " nOrbs=", Method_Mb_OrbBased_nOrbs(1), " nGrid=", Grid_nPoints, " nManyBody=", &
    merge(Coeffs_nCoeffs, TwoRdm_nElements, ciQ), "  total=", tTotal, "  manyBodyPart=", tPart, &
    "  gridPart=", tTotal - tPart

  call Say_Goodbye

end program
