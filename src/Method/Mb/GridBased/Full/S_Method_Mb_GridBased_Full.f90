! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_GridBased_Full) S_Method_Mb_GridBased_Full

  implicit none

  integer(I32), allocatable :: iCompression(:)

  !=============================================================================
  ! local data
  !=============================================================================

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_GridBased_Full_Fabricate
    use M_Utils_Say
    use M_Method

    call Say_Fabricate("method.mb.gridBased.full")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Method_Setup => Setup
    Method_TimeDerivative => TimeDerivative

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Method
    use M_Method_Mb
    use stdlib_sorting

    integer(I32) :: c, bt, i, e, s, size
    integer(I32) :: code(Method_Mb_nBodiesSum)

    call Say_Setup("method.mb.gridBased.full")

    allocate (iCompression(Grid_nPoints**Method_Mb_nBodiesSum))

    size = 0
    do i = 1, Grid_nPoints**Method_Mb_nBodiesSum
      call Decode(code, i)
      do bt = 1, Method_Mb_nBodyTypes
        s = Method_Mb_nBodiesStart(bt)
        e = Method_Mb_nBodiesEnd(bt)
        if (all(code(s:e - 1) < code(s + 1:e))) size = size + 1
        call Sort(code(s:e))
      end do
      call Encode(c, code)
      iCompression(i) = c
    end do

    allocate (Method_state(1:size))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivative(dState, state, time)
    use M_Utils_Constants
    use M_Grid
    use M_Method_Mb
    use M_Method_Mb_GridBased

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64) :: hSlice(Grid_nPoints)
    complex(R64) :: slice(Grid_nPoints)

    complex(R64) :: hSlice2d(Grid_nPoints, Grid_nPoints)
    complex(R64) :: slice2d(Grid_nPoints, Grid_nPoints)

    complex(R64) :: dStateFull(Grid_nPoints**Method_Mb_nBodiesSum)
    complex(R64) :: stateFull(Grid_nPoints**Method_Mb_nBodiesSum)

    integer(I32) :: iBo, pBo, rBo, iBo1, pBo1, rBo1, iBo2, pBo2, rBo2
    integer(I32) :: i1, i2, i, l, m, c, s, k, nG

    nG = Grid_nPoints

    dStateFull(:) = 0.0_R64

    do i = 1, nG**Method_Mb_nBodiesSum
      stateFull(i) = state(iCompression(i))
    end do

    do iBo = 1, Method_Mb_nBodiesSum
      pBo = nG**(iBo - 1)
      rBo = nG**(Method_Mb_nBodiesSum - iBo)
      k = 0
      do l = 1, pBo
        k = k + 1
        c = k
        do m = 1, rBo

          s = c
          do i = 1, nG
            slice(i) = stateFull(c)
            c = c + pBo
          end do

          hSlice(:) = 0.0_R64
          call Method_Mb_GridBased_ApplyKineticOp(hSlice, slice, time, bt_=Method_Mb_bodyTypeOfBody(iBo))
          call Method_Mb_GridBased_ApplyPotentialOp(hSlice, slice, time, bt_=Method_Mb_bodyTypeOfBody(iBo))

          c = s
          do i = 1, nG
            dStateFull(c) = dStateFull(c) + hSlice(i)
            c = c + pBo
          end do

        end do
      end do
    end do

    do iBo2 = 1, Method_Mb_nBodiesSum
      do iBo1 = 1, iBo2 - 1
        pBo2 = nG**(iBo2 - 1)
        pBo1 = nG**(iBo1 - 1)
        rBo2 = nG**(Method_Mb_nBodiesSum - iBo2)
        rBo1 = nG**(Method_Mb_nBodiesSum - iBo1)
        k = 0

        do l = 1, pBo1
          k = k + 1
          c = k

          do m = 1, rBo2

            s = c
            do i2 = 1, nG
              do i1 = 1, nG
                slice2d(i1, i2) = stateFull(c)
                c = c + pBo1
              end do
              c = c + pBo2 - pBo1 * nG
            end do

            hSlice2d(:, :) = 0.0_R64
            call Method_Mb_GridBased_ApplyInteractionOp(hSlice2d, &
                                                        slice2d, &
                                                        Method_Mb_bodyTypeOfBody(iBo1), &
                                                        Method_Mb_bodyTypeOfBody(iBo2), &
                                                        time)

            c = s
            do i2 = 1, nG
              do i1 = 1, nG
                dStateFull(c) = dStateFull(c) + hSlice2d(i1, i2)
                c = c + pBo1
              end do
              c = c + pBo2 - pBo1 * nG
            end do

          end do
        end do
      end do
    end do

    do c = 1, nG**Method_Mb_nBodiesSum
      dState(iCompression(c)) = dStateFull(c)
    end do

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pure subroutine Decode(code, index)
    use M_Grid

    integer(I32), intent(out), contiguous :: code(:)
    integer(I32), intent(in) :: index

    integer(I32) :: ibt, iTmp

    iTmp = index - 1
    do ibt = 1, size(code)

      ! This multibase decomposition extracts the ibt component from the encoded coeff array
      ! element index
      code(ibt) = mod(iTmp, Grid_nPoints) + 1
      iTmp = iTmp / Grid_nPoints

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pure subroutine Encode(index, code)
    use M_Grid

    integer(I32), intent(out) :: index
    integer(I32), intent(in), contiguous  :: code(:)

    integer(I32) :: ibt, iTmp

    iTmp = 1
    index = 0
    do ibt = 1, size(code)
      index = index + (code(ibt) - 1) * iTmp
      iTmp = iTmp * Grid_nPoints
    end do

    index = index + 1

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
