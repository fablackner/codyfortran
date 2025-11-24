! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Centralized kind/type definitions for numeric precision.
!>
!> Exposes named kind parameters for signed integers (I08/I16/I32/I64) and
!> reals/complexes (R32/R64/R128) derived from `iso_fortran_env`.
!> All other modules `use` these to ensure consistent precision choices.
module M_Utils_Types
  use, intrinsic :: iso_fortran_env ! provides various informaton of the environment

  !==========================================================
  ! Named integer constants defining the kind of integer and real
  !==========================================================

  integer, parameter :: I08 = INT8   !selected_Int_kind(2)
  integer, parameter :: I16 = INT16  !selected_Int_kind(4)
  integer, parameter :: I32 = INT32  !selected_Int_kind(9)
  integer, parameter :: I64 = INT64  !selected_Int_kind(19)

  integer, parameter :: R32 = REAL32   ! kind(1.0)   selected_Real_kind(p=6,r=37)
  integer, parameter :: R64 = REAL64   ! kind(1.0_R64) selected_Real_kind(p=14,r=300)
  integer, parameter :: R128 = REAL128 !             selected_Real_kind(p=32, r=4900)

end module
