! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Submodule S_CoeffsInit_Load implements the "Load" CI coefficient
!> initializer that reads a pre-computed coefficient vector from disk.
!>
!> Use Cases
!> ---------
!>   - Restart simulations from a previously saved state
!>   - Initialize with externally computed CI coefficients
!>   - Resume time-dependent propagation from a checkpoint
!>
!> File Format
!> -----------
!> The input file `coeffs.in` must be a binary file containing the complex
!> CI coefficient vector in native Fortran storage order. The file size
!> must exactly match `size(coeffs) * storage_size(coeffs) / 8` bytes.
!>
!> The binary format is produced by `SaveData` from M_Utils_DataStorage,
!> which writes raw complex(R64) arrays without headers.
!>
!> JSON Configuration
!> ------------------
!>   "coeffsInit": { "load": { } }
!>
!> No additional parameters are currently supported; the filename is
!> hardcoded to `coeffs.in` in the working directory.
!>
!> @warning Dimension mismatch between the file and the allocated CI space
!>          will cause undefined behavior or a runtime error from LoadData.
submodule(M_CoeffsInit_Load) S_CoeffsInit_Load

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Registers the load initializer by binding `CoeffsInit_Initialize`
  !> in M_CoeffsInit to the local `Initialize` procedure.
  !>
  !> @note No setup phase is needed; file reading occurs at initialization.
  module subroutine CoeffsInit_Load_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_CoeffsInit

    call Say_Fabricate("coeffsInit.load")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    CoeffsInit_Initialize => Initialize

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Loads CI coefficients from the binary file `coeffs.in`.
  !>
  !> The entire coefficient vector is read in a single I/O operation.
  !> The caller must ensure that `coeffs` is allocated with the correct
  !> dimension matching the stored data.
  !>
  !> @param[out] coeffs  Complex CI coefficient vector to be filled from file.
  subroutine Initialize(coeffs)
    use M_Utils_DataStorage

    complex(R64), intent(out), contiguous :: coeffs(:)

    character(len=*), parameter :: filename = "coeffs.in"

    call LoadData(coeffs, filename, storage_size(coeffs))

  end subroutine

end submodule
