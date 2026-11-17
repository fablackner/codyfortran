!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Check shared by the Ne MCTDHF laser tests with declared orbital m,
! T_Ne3d_11 and T_Ne3d_13 (included after `contains` next to
! Ne3dMctdhxLaser.inc).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Number of nonzero entries of the orbitals in channels other than their
!> declared m.
function OtherMEntries(orbs) result(res)
  use M_Grid_Ylm
  use M_Method_Mb_OrbBased

  complex(R64), intent(in), contiguous :: orbs(:, :)
  integer(I32) :: res

  integer(I32) :: i

  res = 0
  do i = 1, size(orbs, 2)
    res = res + count(Grid_Ylm_mCoord .ne. Method_Mb_OrbBased_mOfOrb(i) .and. abs(orbs(:, i)) > 0.0_R64)
  end do

end function
