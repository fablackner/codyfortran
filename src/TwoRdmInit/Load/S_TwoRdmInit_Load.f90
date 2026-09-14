! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_TwoRdmInit_Load.f90
!> @brief Implementation of file-based 2-RDM initialization.
!>
!> @details
!> Loads 2-RDM tensor blocks from binary files. Each block is stored in
!> a separate file following the naming convention `rdm2BtXX_YY.in`.
!>
!> @par File Discovery
!> The filename is constructed at runtime:
!> @code{.f90}
!> write(filename, '("rdm2Bt",I2.2,"_",I2.2,".in")') bt1, bt2
!> @endcode
!>
!> @par Binary Format
!> Uses `M_Utils_DataStorage.LoadData` which expects unformatted Fortran
!> binary data. The tensor is read directly into the output array with
!> element size determined by `storage_size(rdm2Bt)`.
submodule(M_TwoRdmInit_Load) S_TwoRdmInit_Load

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Bind procedure pointers to file-loading implementation.
  module subroutine TwoRdmInit_Load_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_TwoRdmInit

    call Say_Fabricate("twoRdmInit.load")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Map the local implementation to the main procedure pointer
    TwoRdmInit_InitializeBlock => InitializeBlock

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Load a 2-RDM block from a binary file.
!>
!> @details
!> Constructs the filename from body-type indices and loads the tensor
!> data using the DataStorage utility. The file must exist and contain
!> data compatible with the expected array dimensions.
!>
!> @par Filename Convention
!> @code
!> rdm2Bt{bt1:02d}_{bt2:02d}.in
!> @endcode
!> Examples: rdm2Bt01_01.in, rdm2Bt01_02.in, rdm2Bt02_01.in, rdm2Bt02_02.in
!>
!> @param[out] rdm2Bt  Output tensor (shape must match file contents)
!> @param[in]  bt1     Body type index 1 (used in filename)
!> @param[in]  bt2     Body type index 2 (used in filename)
!>
!> @warning Will crash if file does not exist or has incompatible size.
  subroutine InitializeBlock(rdm2Bt, bt1, bt2)
    use M_Utils_DataStorage

    complex(R64), intent(out), contiguous  :: rdm2Bt(:, :, :, :)
    integer(I32), intent(in)               :: bt1
    integer(I32), intent(in)               :: bt2

    character(len=256) :: filename

    write (filename, '("rdm2Bt",I2.2,"_",I2.2,".in")') bt1, bt2

    call LoadData(rdm2Bt, trim(filename), storage_size(rdm2Bt))

  end subroutine

end submodule
