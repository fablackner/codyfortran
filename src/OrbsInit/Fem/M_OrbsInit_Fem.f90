! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Orbital initialization on the 3D Cartesian spectral-element grid
!> (grid.fem).
!>
!> @details
!> Two-level pattern (as OrbsInit_Linear): the family-level InitializeOrb
!> samples the model's pointwise function
!>   OrbsInit_Fem_InitFunction(x, y, z, ind, bt_)
!> at every interior node and normalizes the result with Grid_InnerProduct.
!> It does not orthonormalize (programs do that after Method_Setup).
!>
!> Available models
!> ----------------
!> - **sto**: linear combinations of Slater-type orbitals with real solid
!>   harmonics on the grid.fem centers (molecular LCAO guesses)
!>
!> @see M_OrbsInit_Fem_Sto, M_Grid_Fem
module M_OrbsInit_Fem
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Check grid.fem, bind OrbsInit_InitializeOrb and dispatch to the model.
    module subroutine OrbsInit_Fem_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Number of orbitals the active model defines; InitializeOrb stops for an
  !> orbital index outside 1..OrbsInit_Fem_nOrbitals (set by the model).
  integer(I32) :: OrbsInit_Fem_nOrbitals = 0

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointwise initializer f(x, y, z; ind, bt_) bound by the model.
  procedure(I_OrbsInit_Fem_InitFunction), pointer :: OrbsInit_Fem_InitFunction => null()
  abstract interface
    !> Real (unnormalized) amplitude of orbital ind at the point (x, y, z).
    !> Called concurrently from the OpenMP node loop of InitializeOrb, so it
    !> must be thread-safe (read-only module data).
    function I_OrbsInit_Fem_InitFunction(x, y, z, ind, bt_) result(res)
      import :: I32, R64
      !> Orbital amplitude at (x, y, z).
      real(R64)                          :: res
      !> Cartesian coordinates of the node.
      real(R64), intent(in)              :: x
      real(R64), intent(in)              :: y
      real(R64), intent(in)              :: z
      !> Orbital index within the body type (1..OrbsInit_Fem_nOrbitals).
      integer(I32), intent(in)           :: ind
      !> Optional body type (absent for method.sb).
      integer(I32), intent(in), optional :: bt_
    end function
  end interface

end module
