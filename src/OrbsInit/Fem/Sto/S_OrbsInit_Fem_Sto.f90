! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for the Slater-type-orbital initializer.
submodule(M_OrbsInit_Fem_Sto) S_OrbsInit_Fem_Sto

  implicit none

  !-----------------------------------------------------------------------------
  ! local data
  !-----------------------------------------------------------------------------

  !> Largest l of all terms (size of the solid harmonic tables).
  integer(I32) :: lmaxTerms = 0

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Read the term arrays, validate them and bind the InitFunction.
  module subroutine OrbsInit_Fem_Sto_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid_Fem
    use M_OrbsInit_Fem

    character(len=200) :: message
    integer(I32) :: nTerms, iTerm, iOrb

    call Say_Fabricate("orbsInit.fem.sto")

    !------------------------------------
    ! read the parallel term arrays
    !------------------------------------

    call RequireKey("orbital")
    call RequireKey("center")
    call RequireKey("n")
    call RequireKey("l")
    call RequireKey("m")
    call RequireKey("zeta")

    OrbsInit_Fem_Sto_orbital = Json_Get("orbsInit.fem.sto.orbital", [integer(I32) ::])
    OrbsInit_Fem_Sto_center = Json_Get("orbsInit.fem.sto.center", [integer(I32) ::])
    OrbsInit_Fem_Sto_n = Json_Get("orbsInit.fem.sto.n", [integer(I32) ::])
    OrbsInit_Fem_Sto_l = Json_Get("orbsInit.fem.sto.l", [integer(I32) ::])
    OrbsInit_Fem_Sto_m = Json_Get("orbsInit.fem.sto.m", [integer(I32) ::])
    OrbsInit_Fem_Sto_zeta = Json_Get("orbsInit.fem.sto.zeta", [real(R64) ::])

    nTerms = size(OrbsInit_Fem_Sto_orbital)
    if (nTerms < 1) error stop "orbsInit.fem.sto needs at least one term"

    OrbsInit_Fem_Sto_coeff = Json_Get("orbsInit.fem.sto.coeff", [(1.0_R64, iTerm=1, nTerms)])

    !------------------------------------
    ! validate
    !------------------------------------

    if (size(OrbsInit_Fem_Sto_center) .ne. nTerms .or. size(OrbsInit_Fem_Sto_n) .ne. nTerms .or. &
        size(OrbsInit_Fem_Sto_l) .ne. nTerms .or. size(OrbsInit_Fem_Sto_m) .ne. nTerms .or. &
        size(OrbsInit_Fem_Sto_zeta) .ne. nTerms .or. size(OrbsInit_Fem_Sto_coeff) .ne. nTerms) then
      error stop "orbsInit.fem.sto: orbital, center, n, l, m, zeta and coeff need one entry per term (equal sizes)"
    end if

    do iTerm = 1, nTerms
      if (OrbsInit_Fem_Sto_orbital(iTerm) < 1) then
        write (message, '(A, I0, A)') "orbsInit.fem.sto: term ", iTerm, " needs orbital >= 1"
        error stop trim(message)
      end if

      if (OrbsInit_Fem_Sto_center(iTerm) < 1 .or. OrbsInit_Fem_Sto_center(iTerm) > Grid_Fem_nCenters) then
        write (message, '(A, I0, A, I0, A)') "orbsInit.fem.sto: term ", iTerm, " needs 1 <= center <= ", &
          Grid_Fem_nCenters, " (the grid.fem centers)"
        error stop trim(message)
      end if

      if (.not. (OrbsInit_Fem_Sto_l(iTerm) >= 0 .and. OrbsInit_Fem_Sto_n(iTerm) > OrbsInit_Fem_Sto_l(iTerm))) then
        write (message, '(A, I0, A)') "orbsInit.fem.sto: term ", iTerm, " needs n > l >= 0"
        error stop trim(message)
      end if

      if (abs(OrbsInit_Fem_Sto_m(iTerm)) > OrbsInit_Fem_Sto_l(iTerm)) then
        write (message, '(A, I0, A)') "orbsInit.fem.sto: term ", iTerm, " needs |m| <= l"
        error stop trim(message)
      end if

      if (.not. (OrbsInit_Fem_Sto_zeta(iTerm) > 0.0_R64)) then
        write (message, '(A, I0, A)') "orbsInit.fem.sto: term ", iTerm, " needs zeta > 0"
        error stop trim(message)
      end if
    end do

    ! Every orbital index up to the largest one needs a term
    do iOrb = 1, maxval(OrbsInit_Fem_Sto_orbital)
      if (.not. any(OrbsInit_Fem_Sto_orbital .eq. iOrb)) then
        write (message, '(A, I0, A)') "orbsInit.fem.sto: orbital ", iOrb, " has no term (indices must cover 1..max)"
        error stop trim(message)
      end if
    end do

    lmaxTerms = maxval(OrbsInit_Fem_Sto_l)

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_Fem_nOrbitals = maxval(OrbsInit_Fem_Sto_orbital)
    OrbsInit_Fem_InitFunction => InitFunction

  contains

    !> Stops unless orbsInit.fem.sto.<key> is given.
    subroutine RequireKey(key)
      character(len=*), intent(in) :: key

      if (.not. Json_GetExistence("orbsInit.fem.sto."//key)) then
        error stop "orbsInit.fem.sto."//key//" is required (one entry per term)"
      end if

    end subroutine

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Amplitude of orbital ind at (x, y, z): the sum of its STO terms
!> coeff rho^(n-l-1) S_lm(r - R) exp(-zeta rho) (no division by rho, regular
!> on the nucleus nodes). The same for every body type; thread-safe.
  function InitFunction(x, y, z, ind, bt_) result(res)
    use M_Utils_UnusedVariables
    use M_Utils_SolidHarmonics, only: SolidHarmonics_FillValues, SolidHarmonics_LmIndex
    use M_Grid_Fem

    real(R64)                          :: res
    real(R64), intent(in)              :: x
    real(R64), intent(in)              :: y
    real(R64), intent(in)              :: z
    integer(I32), intent(in)           :: ind
    integer(I32), intent(in), optional :: bt_

    real(R64) :: solid((lmaxTerms + 1)**2), d(3), rho, radial
    integer(I32) :: iTerm, l, power

    if (.false.) call UnusedVariables_Mark(bt_)

    res = 0.0_R64

    do iTerm = 1, size(OrbsInit_Fem_Sto_orbital)
      if (OrbsInit_Fem_Sto_orbital(iTerm) .ne. ind) cycle

      l = OrbsInit_Fem_Sto_l(iTerm)
      d = [x, y, z] - Grid_Fem_centers(:, OrbsInit_Fem_Sto_center(iTerm))
      rho = norm2(d)

      ! rho^(n-l-1), with rho^0 = 1 also on the nucleus node
      power = OrbsInit_Fem_Sto_n(iTerm) - l - 1
      radial = 1.0_R64
      if (power > 0) radial = rho**power

      call SolidHarmonics_FillValues(solid, d, l)

      res = res + OrbsInit_Fem_Sto_coeff(iTerm) * radial * solid(SolidHarmonics_LmIndex(l, OrbsInit_Fem_Sto_m(iTerm))) * &
            exp(-OrbsInit_Fem_Sto_zeta(iTerm) * rho)
    end do

  end function

end submodule
