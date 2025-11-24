! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Shared metadata and entry points for many-body methods.
!>
!> This module collects common information used by all many-body variants
!> (orbital-, grid-, and GEM-based). It exposes per-body-type counts, index
!> ranges, and statistics needed for building operators and mapping slices of
!> the full state. The `Method_Mb_Fabricate` routine binds method-specific
!> procedure pointers and initializes these arrays from the input configuration.
module M_Method_Mb
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds many-body method procedure pointers at runtime and initializes
    !> shared metadata (body types, counts, index ranges, statistics) from the
    !> input configuration.
    module subroutine Method_Mb_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Number of distinct body types present in the system.
  !> Body types can represent different particle species or spin states.
  !> Must be >= 1.
  integer(I32) :: Method_Mb_nBodyTypes

  !> Classification of particles by body type (length = total number of bodies).
  !> Enables type-dependent potentials/interactions and statistics.
  !> Invariant: 1 <= Method_Mb_bodyTypeOfBody(i) <= Method_Mb_nBodyTypes.
  integer(I32), allocatable :: Method_Mb_bodyTypeOfBody(:)

  !> Quantum statistics per body type ('f' for fermions, 'b' for bosons).
  !> Length = Method_Mb_nBodyTypes. Case-insensitive.
  character, allocatable :: Method_Mb_bodyStatistics(:)

  !> Number of particles per body type.
  !> Method_Mb_nBodies(bt) gives the number of particles of body type bt.
  !> Length = Method_Mb_nBodyTypes. Non-negative.
  integer(I32), allocatable :: Method_Mb_nBodies(:)

  !> 1-based start index of particles per body type in a packed, contiguous layout.
  !> Method_Mb_nBodiesStart(bt) gives the first global body index of type bt.
  integer(I32), allocatable :: Method_Mb_nBodiesStart(:)

  !> 1-based end index of particles per body type in a packed, contiguous layout.
  !> Method_Mb_nBodiesEnd(bt) gives the last global body index of type bt.
  integer(I32), allocatable :: Method_Mb_nBodiesEnd(:)

  !> Total number of particles across all body types.
  !> Invariant: Method_Mb_nBodiesSum = sum(Method_Mb_nBodies).
  integer(I32) :: Method_Mb_nBodiesSum

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
