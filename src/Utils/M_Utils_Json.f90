! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Minimal JSON reading/writing helpers for configuration and results.
!>
!> Provides convenience routines to parse, query, and serialize simple JSON
!> structures used by the app and test harnesses.
module M_Utils_Json
  use M_Utils_Types
  use M_Utils_Characters
  use JSON_MODULE

  implicit none

  type(JSON_FILE)               :: jsonFile

  interface Json_Get

    procedure Json_GetInteger
    procedure Json_GetIntegerArray

    procedure Json_GetReal
    procedure Json_GetRealArray

    procedure Json_GetBool
    procedure Json_GetBoolArray

    procedure Json_GetString
    procedure Json_GetStringArray

  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Json_LoadJsonFile(manualFileName_, relativeToProjectDirQ_)

    character(len=*), intent(in), optional :: manualFileName_
    logical, intent(in), optional          :: relativeToProjectDirQ_

    integer(I32) ::  istat
    character(len=:), allocatable :: jsonFileName
    logical :: fileExists
    logical :: relativeToProjectDirQ

    if (present(relativeToProjectDirQ_)) relativeToProjectDirQ = relativeToProjectDirQ_
    if (.not. present(relativeToProjectDirQ_)) relativeToProjectDirQ = .false.

    if (present(manualFileName_)) then
      jsonFileName = manualFileName_
      if (relativeToProjectDirQ) then
        block
          character(len=255) :: codyProjectDir
          integer :: envLen
          call get_environment_variable("CODY_PROJECT_DIR", value=codyProjectDir, length=envLen, status=istat)
          if (istat .ne. 0 .or. envLen .eq. 0) then
            error stop "CODY_PROJECT_DIR environment variable is not set or empty or too long."
          end if
          jsonFileName = trim(codyProjectDir)//"/"//trim(jsonFileName)
        end block
      end if
    end if

    if (.not. present(manualFileName_)) then
      if (relativeToProjectDirQ) error stop 'relativeToProjectDirQ_ and manualFileName_ are exclusive'

      allocate (character(len=255) :: jsonFileName)
      call get_command_argument(1, jsonFileName, status=istat)
      if (istat .ne. 0) then
        call get_command_argument(0, jsonFileName, status=istat)
        if (istat .ne. 0) error stop 'cannot extract jsonFileName'
        jsonFileName = jsonFileName(:scan(jsonFileName, '.', .true.))//'json'
      end if
    end if

    inquire (file=trim(jsonFileName), exist=fileExists)
    if (.not. fileExists) then
      error stop 'Json_LoadJsonFile: file does not exist -> '//trim(jsonFileName)
    end if

    ! Initialise the JSON_FILE object.
    call jsonFile % initialize()

    ! Load the file.
    call jsonFile % load(filename=trim(jsonFileName))
    if (jsonFile % failed()) error stop 'jsonFileName does not exist or has no proper JSON structure'

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Json_CreateJsonFile(jsonString, manualFileName_)

    character(len=*), intent(in) :: jsonString
    character(len=*), intent(in) :: manualFileName_

    integer :: fd, iostat

    open (newunit=fd, file=manualFileName_, status='replace', action='write', iostat=iostat)
    if (iostat .ne. 0) then
      print *, 'Error opening file.'
      stop
    end if

    ! Write the JSON string to the file
    write (fd, '(a)') trim(jsonString)
    close (fd)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetInteger(name, def, alternativeQ_, path_) result(res)

    integer(I32) :: res
    character(len=*), intent(in)           :: name
    integer(I32), intent(in)               :: def
    logical, intent(out), optional         :: alternativeQ_
    character(len=*), intent(in), optional :: path_

    character(len=:), allocatable          :: path
    logical                                :: foundQ

    ! Read in the data.
    if (present(path_)) path = path_//"."//name !there is a bug in -Wmaybe-uninitialized
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, res, foundQ)

    if (foundQ) then
      write (*, '(1X, A, 3A, *(G0,1X), A)') green, 'in: ', path, ' = ', res, reset
    else

      if (.not. present(alternativeQ_)) then
        res = def
        write (*, '(1X, A, 3A, *(G0,1X), A)') red, 'de: ', path, ' = ', res, reset
      else
        ! when alternativeQ_ is present an alternative definition must be
      end if
    end if

    if (present(alternativeQ_)) alternativeQ_ = .not. foundQ

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetIntegerArray(name, def, alternativeQ_, path_) result(res)

    integer(I32), allocatable              :: res(:)
    character(len=*), intent(in)           :: name
    integer(I32), intent(in), contiguous   :: def(:)
    logical, intent(out), optional         :: alternativeQ_
    character(len=*), intent(in), optional :: path_

    character(len=:), allocatable          :: path
    logical                                :: foundQ

    ! Read in the data.
    if (present(path_)) path = path_//"."//name !there is a bug in -Wmaybe-uninitialized
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, res, foundQ)

    if (foundQ) then
      write (*, '(1X, A, 3A, *(G0,1X), A)') green, 'in: ', path, ' = ', res, reset
    else
      if (.not. present(alternativeQ_)) then
        res = def
        write (*, '(1X, A, 3A, *(G0,1X), A)') red, 'de: ', path, ' = ', res, reset
      else
        ! when alternativeQ_ is present an alternative definition must be
      end if
    end if

    if (present(alternativeQ_)) alternativeQ_ = .not. foundQ

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetReal(name, def, alternativeQ_, path_) result(res)

    real(R64)                              :: res
    character(len=*), intent(in)           :: name
    real(R64), intent(in)                  :: def
    logical, intent(out), optional         :: alternativeQ_
    character(len=*), intent(in), optional :: path_

    character(len=:), allocatable          :: path
    logical                                :: foundQ

    ! Read in the data.
    if (present(path_)) path = path_//"."//name !there is a bug in -Wmaybe-uninitialized
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, res, foundQ)

    if (foundQ) then
      write (*, '(1X, A, 3A, *(G0,1X), A)') green, 'in: ', path, ' = ', res, reset
    else
      if (.not. present(alternativeQ_)) then
        res = def
        write (*, '(1X, A, 3A, *(G0,1X), A)') red, 'de: ', path, ' = ', res, reset
      else
        ! when alternativeQ_ is present an alternative definition must be
      end if
    end if

    if (present(alternativeQ_)) alternativeQ_ = .not. foundQ

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetRealArray(name, def, alternativeQ_, path_) result(res)

    real(R64), allocatable                 :: res(:)
    character(len=*), intent(in)           :: name
    real(R64), intent(in), contiguous      :: def(:)
    logical, intent(out), optional         :: alternativeQ_
    character(len=*), intent(in), optional :: path_

    character(len=:), allocatable          :: path
    logical                                :: foundQ

    ! Read in the data.
    if (present(path_)) path = path_//"."//name !there is a bug in -Wmaybe-uninitialized
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, res, foundQ)

    if (foundQ) then
      write (*, '(1X, A, 3A, *(G0,1X), A)') green, 'in: ', path, ' = ', res, reset
    else
      if (.not. present(alternativeQ_)) then
        res = def
        write (*, '(1X, A, 3A, *(G0,1X), A)') red, 'de: ', path, ' = ', res, reset
      else
        ! when alternativeQ_ is present an alternative definition must be
      end if
    end if

    if (present(alternativeQ_)) alternativeQ_ = .not. foundQ

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetBool(name, def, alternativeQ_, path_) result(res)

    logical                                :: res
    character(len=*), intent(in)           :: name
    logical, intent(in)                    :: def
    logical, intent(out), optional         :: alternativeQ_
    character(len=*), intent(in), optional :: path_

    character(len=:), allocatable          :: path
    logical                                :: foundQ

    ! Read in the data.
    if (present(path_)) path = path_//"."//name !there is a bug in -Wmaybe-uninitialized
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, res, foundQ)

    if (foundQ) then
      write (*, '(1X, A, 3A, *(G0,1X), A)') green, 'in: ', path, ' = ', res, reset
    else
      if (.not. present(alternativeQ_)) then
        res = def
        write (*, '(1X, A, 3A, *(G0,1X), A)') red, 'de: ', path, ' = ', res, reset
      else
        ! when alternativeQ_ is present an alternative definition must be
      end if
    end if

    if (present(alternativeQ_)) alternativeQ_ = .not. foundQ

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetBoolArray(name, def, alternativeQ_, path_) result(res)

    logical, allocatable                   :: res(:)
    character(len=*), intent(in)           :: name
    logical, intent(in), contiguous        :: def(:)
    logical, intent(out), optional         :: alternativeQ_
    character(len=*), intent(in), optional :: path_

    character(len=:), allocatable          :: path
    logical                                :: foundQ

    ! Read in the data.
    if (present(path_)) path = path_//"."//name !there is a bug in -Wmaybe-uninitialized
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, res, foundQ)

    ! Read in the data.
    if (foundQ) then
      write (*, '(1X, A, 3A, *(G0,1X), A)') green, 'in: ', path, ' = ', res, reset
    else
      if (.not. present(alternativeQ_)) then
        res = def
        write (*, '(1X, A, 3A, *(G0,1X), A)') red, 'de: ', path, ' = ', res, reset
      else
        ! when alternativeQ_ is present an alternative definition must be
      end if
    end if

    if (present(alternativeQ_)) alternativeQ_ = .not. foundQ

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetString(name, def, alternativeQ_, path_) result(res)

    character(len=:), allocatable          :: res
    character(len=*), intent(in)           :: name
    character(len=*), intent(in)           :: def
    logical, intent(out), optional         :: alternativeQ_
    character(len=*), intent(in), optional :: path_

    character(len=:), allocatable          :: path
    logical                                :: foundQ

    ! Read in the data.
    if (present(path_)) path = path_//"."//name !there is a bug in -Wmaybe-uninitialized
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, res, foundQ)

    if (foundQ) then
      write (*, '(1X, A, 3A, *(G0,1X), A)') green, 'in: ', path, ' = ', res, reset
    else
      if (.not. present(alternativeQ_)) then
        res = def
        write (*, '(1X, A, 3A, *(G0,1X), A)') red, 'de: ', path, ' = ', res, reset
      else
        ! when alternativeQ_ is present an alternative definition must be
      end if
    end if

    if (present(alternativeQ_)) alternativeQ_ = .not. foundQ

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetStringArray(name, def, alternativeQ_, path_) result(res)

    character(len=:), allocatable          :: res(:)
    character(len=*), intent(in)           :: name
    character(len=*), intent(in), contiguous :: def(:)
    logical, intent(out), optional         :: alternativeQ_
    character(len=*), intent(in), optional :: path_

    character(len=:), allocatable          :: path
    integer(I32), allocatable              :: stringLengths(:)
    logical                                :: foundQ

    if (present(path_)) path = path_//"."//name
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, res, stringLengths, foundQ)

    if (foundQ) then
      write (*, '(1X, A, 3A, *(G0,1X), A)') green, 'in: ', path, ' = ', res, reset
    else
      if (.not. present(alternativeQ_)) then
        res = def
        write (*, '(1X, A, 3A, *(G0,1X), A)') red, 'de: ', path, ' = ', res, reset
      else
        allocate (character(len=0) :: res(0))
      end if
    end if

    if (present(alternativeQ_)) alternativeQ_ = .not. foundQ

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetExistence(name, path_) result(res)

    logical                                :: res
    character(len=*), intent(in)           :: name
    character(len=*), intent(in), optional :: path_

    type(JSON_CORE)                        :: jsonCore
    type(JSON_VALUE), POINTER              :: object
    character(len=:), allocatable          :: path
    logical                                :: foundQ

    call jsonFile % Get_Core(jsonCore)

    ! Read in the data.
    if (present(path_)) path = path_//"."//name !there is a bug in -Wmaybe-uninitialized
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, object, foundQ)

    res = foundQ

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetNumChildren(name, path_) result(res)

    integer(I32)                           :: res
    character(len=*), intent(in)           :: name
    character(len=*), intent(in), optional :: path_

    type(JSON_CORE)                        :: jsonCore
    type(JSON_VALUE), POINTER              :: object
    character(len=:), allocatable          :: path
    logical                                :: foundQ

    call jsonFile % Get_Core(jsonCore)

    ! Read in the data.
    if (present(path_)) path = path_//"."//name !there is a bug in -Wmaybe-uninitialized
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, object, foundQ)
    !if (.not. foundQ) error stop 'res: object not found'

    if (.not. foundQ) then
      res = -1
    else
      res = jsonCore % count(object)
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetChildName(name, index, path_) result(res)

    character(len=:), allocatable          :: res
    character(len=*), intent(in)           :: name
    integer(I32), intent(in)               :: index
    character(len=*), intent(in), optional :: path_

    type(JSON_CORE)                        :: jsonCore
    type(JSON_VALUE), POINTER              :: object
    type(JSON_VALUE), POINTER              :: child
    character(len=:), allocatable          :: path
    logical                                :: foundQ

    call jsonFile % Get_Core(jsonCore)

    ! Read in the data.
    if (present(path_)) path = path_//"."//name !there is a bug in -Wmaybe-uninitialized
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, object, foundQ)
    if (.not. foundQ) error stop 'GetOwnPropertyNames: object not found'

    call jsonCore % Get_Child(object, index, child, foundQ)
    if (.not. foundQ) error stop 'GetOwnPropertyNames: child not found'

    call jsonCore % info(child, name=res)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Json_GetChildType(name, index, path_) result(res)

    character(len=:), allocatable          :: res
    character(len=*), intent(in)           :: name
    integer(I32), intent(in)               :: index
    character(len=*), intent(in), optional :: path_

    integer(I32)                           :: iType
    type(JSON_CORE)                        :: jsonCore
    type(JSON_VALUE), POINTER              :: object
    type(JSON_VALUE), POINTER              :: child
    character(len=:), allocatable          :: path
    logical                                :: foundQ

    call jsonFile % Get_Core(jsonCore)

    ! Read in the data.
    if (present(path_)) path = path_//"."//name !there is a bug in -Wmaybe-uninitialized
    if (.not. present(path_)) path = name

    call jsonFile % Get(path, object, foundQ)
    if (.not. foundQ) error stop 'GetOwnPropertyNames: object not found'

    call jsonCore % Get_Child(object, index, child, foundQ)
    if (.not. foundQ) error stop 'GetOwnPropertyNames: child not found'

    call jsonCore % info(child, var_type=iType)

    select case (iType)

    case (0)
      res = "unknown"
    case (1)
      res = "null"
    case (2)
      res = "object"
    case (3)
      res = "array"
    case (4)
      res = "logical"
    case (5)
      res = "integer"
    case (6)
      res = "real"
    case (7)
      res = "string"
    case (8)
      res = "double" ! equivalent to real

    end select

  end function

end module
