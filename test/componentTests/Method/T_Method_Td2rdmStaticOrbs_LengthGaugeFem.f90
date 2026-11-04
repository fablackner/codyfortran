!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> TD-2RDM with static orbitals (method.mb.orbBased.td2rdmStaticOrbs) and the
!> length-gauge laser of grid.fem (sysGauge.fem.length): the one-body
!> Hamiltonian h1 of the static orbitals follows the field, h2 stays static.
!> Two electrons at the proton of T_Method_Tdci_LengthGaugeFem (same mesh,
!> same pulse along z) with the Coulomb interaction: restricted orbitals 1s
!> (zeta 1) and 2p_z (zeta 0.5) of orbsInit.fem.sto, singlet 2-RDM of the 1s^2
!> determinant (twoRdmInit.excited), Valdemoro reconstruction.
!>
!> td2rdmStaticOrbs fills h1 and h2 at setup (t = 0, before the pulse) and
!> refills h1 in Method_TimeDerivative when a term of it depends on time, now
!> including a configured gauge that is not SysGauge_timeIndependentQ (as
!> TDCI; before, h1 stayed frozen at t = 0 for every SysGauge laser).
!>
!> Checks:
!> - Method_TimeDerivative before the pulse (t = 0): the refilled h1 equals
!>   that of the setup (relative 1e-14)
!> - Method_TimeDerivative at the envelope peak (E = F0 = 0.05):
!>   h1(t) - h1(0) = E(t) D in both spin blocks (restricted: the spatial block
!>   twice) and zero between them, with the dipole matrix D_ij = <phi_i|z|phi_j>
!>   computed here (relative 1e-12)
!> - h2 unchanged by both time derivatives (static orbitals and interaction)
!> - the coupling is not trivial: |D(1s, 2p_z)| > 0.1
program T_Method_Td2rdmStaticOrbs_LengthGaugeFem
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use M_SysKinetic
  use M_SysPotential
  use M_SysGauge
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_TwoRdm
  use M_TwoRdmInit
  use M_Reconstruction
  use M_Method
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  complex(R64), allocatable :: h1Start(:, :), h1Before(:, :), h1Peak(:, :), h2Start(:, :, :, :), dipole(:, :)
  complex(R64), allocatable :: dipoleSpin(:, :), zOrb(:), dState(:)
  real(R64) :: fieldStrength, omega, nCycles, cep, tStart, duration, tPeak, fieldPeak, devBefore, devPeak, devH2
  integer(I32) :: nOS, i, j, ix, iy, iz, node

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Method/T_Method_Td2rdmStaticOrbs_LengthGaugeFem.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysGauge_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call TwoRdm_Fabricate
  call TwoRdmInit_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call Reconstruction_Fabricate

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysGauge_Setup
  call SysInteraction_Setup
  call TwoRdm_Setup
  call TwoRdmInit_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call Method_Setup
  call Reconstruction_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  fieldStrength = Json_Get("sysGauge.fem.length.fieldStrength", 0.0_R64)
  omega = Json_Get("sysGauge.fem.length.omega", 1.0_R64)
  nCycles = Json_Get("sysGauge.fem.length.nCycles", 1.0_R64)
  cep = Json_Get("sysGauge.fem.length.cep", 0.0_R64)
  tStart = Json_Get("sysGauge.fem.length.tStart", 0.0_R64)
  duration = 2.0_R64 * PI * nCycles / omega
  tPeak = tStart + 0.5_R64 * duration
  fieldPeak = fieldStrength * sin(PI * (tPeak - tStart) / duration)**2 * &
              cos(omega * (tPeak - tStart - 0.5_R64 * duration) + cep)

  ! h1 and h2 of the setup (t = 0), h1 after a time derivative before the
  ! pulse and after one at the envelope peak
  h1Start = Method_Mb_OrbBased_h1
  h2Start = Method_Mb_OrbBased_h2
  allocate (dState, mold=Method_state)
  call Method_TimeDerivative(dState, Method_state, 0.0_R64)
  h1Before = Method_Mb_OrbBased_h1
  call Method_TimeDerivative(dState, Method_state, tPeak)
  h1Peak = Method_Mb_OrbBased_h1

  ! Dipole matrix of the static (spatial) orbitals, in both spin blocks
  nOS = Orbs_nOrbsInState
  allocate (dipole(nOS, nOS), zOrb(Grid_nPoints))
  do j = 1, nOS
    do iz = 1, Grid_Fem_nz
      do iy = 1, Grid_Fem_ny
        do ix = 1, Grid_Fem_nx
          node = Grid_Fem_Index(ix, iy, iz)
          zOrb(node) = Grid_Fem_axes(3) % points(iz) * Orbs_orbs(node, j)
        end do
      end do
    end do
    do i = 1, nOS
      dipole(i, j) = Grid_InnerProduct(Orbs_orbs(:, i), zOrb)
    end do
  end do
  allocate (dipoleSpin(2 * nOS, 2 * nOS))
  dipoleSpin(:, :) = (0.0_R64, 0.0_R64)
  dipoleSpin(:nOS, :nOS) = dipole(:, :)
  dipoleSpin(nOS + 1:, nOS + 1:) = dipole(:, :)

  devBefore = maxval(abs(h1Before - h1Start)) / maxval(abs(h1Start))
  devPeak = maxval(abs(h1Peak - h1Start - fieldPeak * dipoleSpin)) / (abs(fieldPeak) * maxval(abs(dipole)))
  devH2 = maxval(abs(Method_Mb_OrbBased_h2 - h2Start)) / maxval(abs(h2Start))

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(1X, A, I0, A, I0)') "mesh N = ", Grid_nPoints, ", spatial orbitals ", nOS
  print *, "envelope peak t, E(t):                          ", tPeak, fieldPeak
  print *, "dipole matrix <phi_i|z|phi_j>:                  ", real(dipole)
  print *, "h1(0) (spin orbitals):                          ", real(h1Start)
  print *, "h1(t_peak):                                     ", real(h1Peak)
  print *, "max |h1(t = 0) - h1(setup)| / |h1|:             ", devBefore
  print *, "max |h1(t) - h1(0) - E(t) D| / |E D|:           ", devPeak
  print *, "max |h2 - h2(setup)| / |h2|:                    ", devH2

  call check(error, abs(dipole(1, 2)) > 0.1_R64)
  if (allocated(error)) error stop "T_Method_Td2rdmStaticOrbs_LengthGaugeFem failure: 1s-2p_z dipole coupling too small"

  call check(error, devBefore < 1e-14_R64)
  if (allocated(error)) error stop "T_Method_Td2rdmStaticOrbs_LengthGaugeFem failure: h1 before the pulse is not h1(0)"

  call check(error, devPeak < 1e-12_R64)
  if (allocated(error)) error stop "T_Method_Td2rdmStaticOrbs_LengthGaugeFem failure: h1 at the peak is not h1(0) + E(t) D"

  call check(error, devH2 < 1e-15_R64)
  if (allocated(error)) error stop "T_Method_Td2rdmStaticOrbs_LengthGaugeFem failure: h2 changed in the time derivative"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
