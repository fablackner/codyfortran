!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Check shared by the Ne MCTDHF laser tests T_Ne3d_07 .. T_Ne3d_12 (included
! after `contains` next to Ne3dMctdhxLaser.inc).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Norm of a packed state in the metric of integrator 1.
function MetricNorm(state) result(res)
  use M_IntegratorList

  complex(R64), intent(in), contiguous :: state(:)
  real(R64) :: res

  res = sqrt(sum(integratorList(1) % e % metricWeights(:) * abs(state(:))**2))

end function
