!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> One- and two-electron integrals of a fixed orbital basis on grid.fem: the
!> Hamiltonian of a finite-basis (static-orbital) calculation, here the 8
!> occupied RHF orbitals of ../../groundstate/scf plus the STOs of
!> P_IntegralsFem.json.
!>
!> The functions come from orbsInit (restricted, nOrbs = number of basis
!> functions; orbsInit.load with nOrbsIn of HF orbitals plus orbsInit.fem.sto
!> for the rest). The integrals are those of the raw, non-orthogonal
!> functions, with their overlap matrix (orthogonalization in the finite
!> basis, make_inputs.py of ../../dynamics/td2rdm). Method td2rdmStaticOrbs
!> fills h1 and h2 (n(n+1)/2 Poisson solves). Written to program.outputFile
!> (default ints.bin, stream access, little endian):
!>   int32 n, real64 nuclear repulsion,
!>   complex128 h1(n, n), h2(n, n, n, n), x(n, n), y(n, n), z(n, n), S(n, n)
!> column-major, spatial: h1(i, j) = <i|T + V|j>, h2(i, j, k, l) =
!> <ij|kl> / 2 (the convention of Method_Mb_OrbBased_FillH2, up-down block),
!> x, y, z the dipole matrices <i|x|j> (coordinates of grid.fem).
program P_IntegralsFem
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Fem
  use M_SysKinetic
  use M_SysPotential
  use M_SysPotential_Fem_Coulomb, only: SysPotential_Fem_Coulomb_nuclearRepulsion
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_TwoRdm
  use M_TwoRdmInit
  use M_Method
  use M_Method_Mb_OrbBased

  implicit none

  complex(R64), allocatable :: dip(:, :, :), rOrbs(:, :), mat(:, :)
  real(R64), allocatable :: coord(:)
  real(R64) :: direction(3)
  character(len=:), allocatable :: outputFile
  integer(I32) :: n, a, i, io

  call Say_Hello
  call Json_LoadJsonFile()

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call TwoRdm_Fabricate
  call TwoRdmInit_Fabricate
  outputFile = Json_Get("program.outputFile", "ints.bin")

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call TwoRdm_Setup
  call TwoRdmInit_Setup
  call Method_Setup

  n = Orbs_nOrbsInState
  allocate (dip(n, n, 3))
  allocate (rOrbs, mold=Orbs_orbs)
  do a = 1, 3
    direction = 0.0_R64
    direction(a) = 1.0_R64
    call Grid_Fem_FillNodeCoordinate(coord, direction)
    do i = 1, n
      rOrbs(:, i) = coord(:) * Orbs_orbs(:, i)
    end do
    call Method_Mb_OrbBased_FillOneBodyMatrix(mat, Orbs_orbs, rOrbs)
    dip(:, :, a) = mat(1:n, 1:n)
  end do

  open (newunit=io, file=outputFile, access="stream", form="unformatted", status="replace")
  write (io) n, SysPotential_Fem_Coulomb_nuclearRepulsion
  write (io) Method_Mb_OrbBased_h1(1:n, 1:n)
  write (io) Method_Mb_OrbBased_h2(1:n, n + 1:2 * n, 1:n, n + 1:2 * n)
  write (io) dip
  call Method_Mb_OrbBased_FillOneBodyMatrix(mat, Orbs_orbs, Orbs_orbs)
  write (io) mat(1:n, 1:n)
  close (io)

  write (*, '(A,I0,A,A)') " integrals of ", n, " orbitals written to ", outputFile
  write (*, '(A,ES12.4)') " max |h2(up-up) - h2(up-down)|: ", &
    maxval(abs(Method_Mb_OrbBased_h2(1:n, 1:n, 1:n, 1:n) - Method_Mb_OrbBased_h2(1:n, n + 1:, 1:n, n + 1:)))

  call Say_Goodbye

end program
