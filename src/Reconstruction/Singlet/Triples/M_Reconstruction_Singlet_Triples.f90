! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Connected triples of the singlet closure (`reconstruction.singlet.triples`).
!>
!> @details
!> An extension of the singlet family: it adds to the closure of the scheme
!> the block of the 3-cumulant that annihilates three occupied and creates
!> three unoccupied orbitals (and its adjoint), which the Nakatsuji-Yasuda
!> products miss. The block follows from the 2-cumulant, the interaction and
!> the orbital energies. A variant binds
!> Reconstruction_Singlet_Triples_FillAmplitudes (amplitudes from the source);
!> `moellerPlesset` is the only one.
!>
!> Layout of the contract: the up-up-down block in the semicanonical
!> orbitals, x(i,j,k,a,b,c) for <a+_a,up a+_b,up a+_c,down a_k,down a_j,up
!> a_i,up> (ket-first), i, j, k occupied, a, b, c unoccupied orbitals. Details:
!> src/Reconstruction/Singlet/Triples/A_Reconstruction_Singlet_Triples.md.
module M_Reconstruction_Singlet_Triples
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Reads reconstruction.singlet.triples, binds the variant, wraps the
    !> closure bound by the singlet family and binds
    !> Reconstruction_SetHamiltonian.
    module subroutine Reconstruction_Singlet_Triples_Fabricate
    end subroutine

    !> The triples block of the 3-cumulant and its adjoint (up-up-down,
    !> ket-first, orbitals of the 2-RDM) for the Hamiltonian of the last
    !> Reconstruction_SetHamiltonian: what the extension adds to the closure.
    module subroutine Reconstruction_Singlet_Triples_FillCumul3(cumul3, rdm1a, rdm2ab)
      complex(R64), intent(out), allocatable :: cumul3(:, :, :, :, :, :)
      complex(R64), intent(in), contiguous   :: rdm1a(:, :)
      complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> The amplitudes of the selected variant (nOcc^3 x nVirt^3) from the source
  !> of the same shape and the orbital energies of the occupied and the
  !> unoccupied semicanonical orbitals.
  procedure(I_Reconstruction_Singlet_Triples_FillAmplitudes), pointer :: &
    Reconstruction_Singlet_Triples_FillAmplitudes => null()
  abstract interface
    subroutine I_Reconstruction_Singlet_Triples_FillAmplitudes(amplitudes, source, epsOcc, epsVirt)
      import :: R64
      complex(R64), intent(out), contiguous :: amplitudes(:, :, :, :, :, :)
      complex(R64), intent(in), contiguous  :: source(:, :, :, :, :, :)
      real(R64), intent(in), contiguous     :: epsOcc(:)
      real(R64), intent(in), contiguous     :: epsVirt(:)
    end subroutine
  end interface

end module
