! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Fem_Laplacian.f90
!> @brief Implementation of the spectral-element kinetic operator and its
!> exact shifted inverse.
submodule(M_SysKinetic_Fem_Laplacian) S_SysKinetic_Fem_Laplacian

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Read and validate the masses, set the flags, bind the operators.
  module subroutine SysKinetic_Fem_Laplacian_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic

    call Say_Fabricate("sysKinetic.fem.laplacian")

    !------------------------------------
    ! read and validate the masses
    !------------------------------------

    SysKinetic_Fem_Laplacian_bodyMass = Json_Get("sysKinetic.fem.laplacian.bodyMass", [1.0_R64])

    if (size(SysKinetic_Fem_Laplacian_bodyMass) < 1) then
      error stop "sysKinetic.fem.laplacian.bodyMass needs at least one entry"
    end if
    if (any(.not. (SysKinetic_Fem_Laplacian_bodyMass > 0.0_R64))) then
      error stop "sysKinetic.fem.laplacian.bodyMass: every mass must be > 0"
    end if

    !------------------------------------
    ! flags and procedure pointers
    !------------------------------------

    ! One mass serves every body type; several masses must cover every body
    ! type a caller requests (MassOfBodyType stops otherwise), and T is body-
    ! type independent only if they are all equal
    SysKinetic_timeIndependentQ = .true.
    SysKinetic_bodyTypeIndependentQ = all(abs(SysKinetic_Fem_Laplacian_bodyMass - &
                                              SysKinetic_Fem_Laplacian_bodyMass(1)) <= 0.0_R64)

    SysKinetic_MultiplyWithKineticOp => MultiplyWithKineticOp
    SysKinetic_SolveShifted => SolveShifted

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief dOrb = T orb = -lap orb / (2 m_bt) (nodal Laplacian, Dirichlet box).
  subroutine MultiplyWithKineticOp(dOrb, orb, time, bt_)
    use M_Utils_UnusedVariables
    use M_Grid_Fem

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous  :: orb(:)
    real(R64), intent(in)                 :: time
    integer(I32), intent(in), optional    :: bt_

    if (.false.) call UnusedVariables_Mark(time)

    call Grid_Fem_ApplyLaplacian(dOrb, orb)
    call ScaleField(dOrb, -0.5_R64 / MassOfBodyType(bt_=bt_))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Exact solution of (T + shift) dOrb = orb for a real or complex
  !> shift:
  !>   (-lap/(2m) + shift) dOrb = orb  <=>  dOrb = 2m (-lap + 2m shift)^-1 orb,
  !> by fast diagonalization (real mode products into and out of the
  !> eigenbasis of -lap; a complex shift mixes the real and imaginary
  !> coefficients only on the diagonal). A real shift must exceed -(smallest
  !> eigenvalue of T).
  subroutine SolveShifted(dOrb, orb, shift, bt_)
    use M_Grid_Fem

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous  :: orb(:)
    complex(R64), intent(in)              :: shift
    integer(I32), intent(in), optional    :: bt_

    real(R64) :: twoMass

    twoMass = 2.0_R64 * MassOfBodyType(bt_=bt_)

    call Grid_Fem_SolveShiftedLaplacian(dOrb, orb, twoMass * shift)
    call ScaleField(dOrb, twoMass)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Mass of body type bt_ (bodyMass(1) if bt_ is absent): a single
  !> bodyMass entry is the mass of every body type; with several entries a
  !> body type outside 1..size(bodyMass) stops (bt_ < 1 always stops).
  function MassOfBodyType(bt_) result(mass)

    integer(I32), intent(in), optional :: bt_
    real(R64) :: mass

    character(len=160) :: message
    integer(I32) :: bt, nMasses

    bt = 1
    if (present(bt_)) bt = bt_
    nMasses = int(size(SysKinetic_Fem_Laplacian_bodyMass), I32)

    if (bt < 1 .or. (nMasses > 1 .and. bt > nMasses)) then
      write (message, '(A, I0, A, I0, A)') "sysKinetic.fem.laplacian: no mass for body type ", bt, &
        " (bodyMass has ", nMasses, " entries; give one mass for all body types or one per body type)"
      error stop trim(message)
    end if

    mass = SysKinetic_Fem_Laplacian_bodyMass(min(bt, nMasses))

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief f = factor f (OpenMP over the nodes outside parallel regions).
  subroutine ScaleField(f, factor)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(inout), contiguous :: f(:)
    real(R64), intent(in)                   :: factor

    integer(I32) :: i

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(f)
      f(i) = factor * f(i)
    end do
    !$omp end parallel do

  end subroutine

end submodule
