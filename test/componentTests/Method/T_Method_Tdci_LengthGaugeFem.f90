!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> TDCI with the length-gauge laser of grid.fem (sysGauge.fem.length): the
!> one-body Hamiltonian h1 of the static orbitals follows the field. H atom on
!> the mesh of T_Method_Sb_LengthGaugeFem with the orbitals 1s (zeta 1) and
!> 2p_z (zeta 0.5) of orbsInit.fem.sto, one electron, the same pulse along z.
!>
!> TDCI fills h1 at setup (t = 0, before the pulse) and refills it in
!> Method_TimeDerivative when a term of it depends on time, now including a
!> configured gauge that is not SysGauge_timeIndependentQ (before, h1 stayed
!> frozen at t = 0 for every SysGauge laser, also the Ylm velocity gauge).
!>
!> Checks:
!> - after Method_TimeDerivative at the envelope peak (E = F0 = 0.05):
!>   h1(t) - h1(0) = E(t) D with the dipole matrix D_ij = <phi_i|z|phi_j>
!>   computed here (relative 1e-12)
!> - the coupling is not trivial: |D(1s, 2p_z)| > 0.1
program T_Method_Tdci_LengthGaugeFem
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Grid
  use M_Grid_Fem
  use M_SysKinetic
  use M_SysPotential
  use M_SysGauge
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  complex(R64), allocatable :: h1Start(:, :), h1Peak(:, :), dipole(:, :), zOrb(:), dState(:)
  real(R64) :: fieldStrength, omega, nCycles, cep, tStart, duration, tPeak, fieldPeak, dev
  integer(I32) :: nOS, i, j, ix, iy, iz, node

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Method/T_Method_Tdci_LengthGaugeFem.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysGauge_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysGauge_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup

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

  ! h1 of the setup (t = 0, before the pulse) and after a time derivative at
  ! the envelope peak
  h1Start = Method_Mb_OrbBased_h1
  allocate (dState, mold=Method_state)
  call Method_TimeDerivative(dState, Method_state, tPeak)
  h1Peak = Method_Mb_OrbBased_h1

  ! Dipole matrix of the static orbitals
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

  dev = maxval(abs(h1Peak - h1Start - fieldPeak * dipole)) / (abs(fieldPeak) * maxval(abs(dipole)))

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(1X, A, I0)') "mesh N = ", Grid_nPoints
  print *, "envelope peak t, E(t):                     ", tPeak, fieldPeak
  print *, "dipole matrix <phi_i|z|phi_j>:             ", real(dipole)
  print *, "h1(0):                                     ", real(h1Start)
  print *, "h1(t_peak):                                ", real(h1Peak)
  print *, "max |h1(t) - h1(0) - E(t) D| / |E D|:      ", dev

  call check(error, abs(dipole(1, 2)) > 0.1_R64)
  if (allocated(error)) error stop "T_Method_Tdci_LengthGaugeFem failure: 1s-2p_z dipole coupling too small to be tested"

  call check(error, dev < 1e-12_R64)
  if (allocated(error)) error stop "T_Method_Tdci_LengthGaugeFem failure: h1 at the envelope peak is not h1(0) + E(t) D"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
