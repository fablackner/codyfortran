program T_Reconstruction_SingletTr0_Ny
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use M_Utils_Types
  use M_Utils_Json
  use M_Utils_RdmEnforceTraceM2ab
  use M_Reconstruction
  use M_Method_Mb
  use M_Method_Mb_OrbBased
  use M_TwoRdm
  use M_TwoRdm_Fermionic_SingletTr0
  use testdrive, only: check, error_type
  implicit none

  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName

  integer(I32) :: nOa, nOb
  integer(I32) :: i1, i2, j1, j2
  complex(R64), allocatable :: rdm2abTr0(:, :, :, :)
  complex(R64), allocatable :: rdm1a(:, :)
  complex(R64), allocatable :: twoRdmState(:)
  complex(R64), allocatable :: dTwoRdmState(:)
  complex(R64), allocatable :: h2(:, :, :, :)

  nOa = 3
  nOb = 3

  Method_Mb_nBodyTypes = 2
  allocate (Method_Mb_nBodies(2))
  Method_Mb_nBodies = [1_I32, 1_I32]
  Method_Mb_nBodiesSum = sum(Method_Mb_nBodies)
  allocate (character(len=1) :: Method_Mb_bodyStatistics(2))
  Method_Mb_bodyStatistics = ['f', 'f']

  allocate (Method_Mb_OrbBased_nOrbs(2))
  Method_Mb_OrbBased_nOrbs = [nOa, nOb]
  Method_Mb_OrbBased_nOrbsSum = nOa + nOb
  allocate (Method_Mb_OrbBased_nOrbsStart(2))
  allocate (Method_Mb_OrbBased_nOrbsEnd(2))
  Method_Mb_OrbBased_nOrbsStart = [1_I32, nOa + 1_I32]
  Method_Mb_OrbBased_nOrbsEnd = [nOa, nOa + nOb]

  call TwoRdm_Fermionic_SingletTr0_Fabricate

  allocate (rdm1a(nOa, nOa))
  rdm1a(:, :) = 0.0_R64
  do i1 = 1, nOa
    rdm1a(i1, i1) = 0.5_R64
  end do

  allocate (rdm2abTr0(nOa, nOb, nOa, nOb))
  do j2 = 1, nOb
    do j1 = 1, nOa
      do i2 = 1, nOb
        do i1 = 1, nOa
          rdm2abTr0(i1, i2, j1, j2) = cmplx(0.01_R64 * (i1 + j1), 0.005_R64 * (i2 - j2), R64)
        end do
      end do
    end do
  end do
  call RdmEnforceTraceM2ab_Apply(rdm2abTr0)

  allocate (twoRdmState(TwoRdm_nElements))
  call TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdms(twoRdmState, rdm2abTr0_=rdm2abTr0, rdm1a_=rdm1a)

  allocate (dTwoRdmState(size(twoRdmState)))
  dTwoRdmState = 0.0_R64

  allocate (h2(nOa + nOb, nOa + nOb, nOa + nOb, nOa + nOb))
  h2(:, :, :, :) = 0.0_R64
  do i1 = 1, nOa + nOb
    h2(i1, i1, i1, i1) = 0.2_R64
  end do

  jsonFileName = "test/componentTests/Reconstruction/SingletTr0/T_Reconstruction_SingletTr0_Ny.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Reconstruction_Fabricate
  call Reconstruction_AddCollisionOp(dTwoRdmState, twoRdmState, h2)

  call check(error, all(.not. ieee_is_nan(real(dTwoRdmState))))
  if (.not. allocated(error)) call check(error, all(.not. ieee_is_nan(aimag(dTwoRdmState))))
  if (allocated(error)) error stop "T_Reconstruction_SingletTr0_Ny failed"

end program T_Reconstruction_SingletTr0_Ny
