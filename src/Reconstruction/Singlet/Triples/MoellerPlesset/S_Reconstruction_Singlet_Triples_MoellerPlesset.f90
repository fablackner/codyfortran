! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> amplitudes = source / D, D = eps_i + eps_j + eps_k - eps_a - eps_b - eps_c:
!> the second-order connected triples (the triples of CCSD(T)) with the
!> 2-cumulant as doubles.
submodule(M_Reconstruction_Singlet_Triples_MoellerPlesset) S_Reconstruction_Singlet_Triples_MoellerPlesset

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Singlet_Triples_MoellerPlesset_Fabricate
    use M_Utils_Say
    use M_Reconstruction_Singlet_Triples

    call Say_Fabricate("reconstruction.singlet.triples.moellerPlesset")

    Reconstruction_Singlet_Triples_FillAmplitudes => FillAmplitudes

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillAmplitudes(amplitudes, source, epsOcc, epsVirt)
    complex(R64), intent(out), contiguous :: amplitudes(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous  :: source(:, :, :, :, :, :)
    real(R64), intent(in), contiguous     :: epsOcc(:)
    real(R64), intent(in), contiguous     :: epsVirt(:)

    integer(I32) :: i, j, k, a, b, c

    !$omp parallel do collapse(2) private(i, j, k, a, b, c)
    do c = 1, size(epsVirt)
      do b = 1, size(epsVirt)
        do a = 1, size(epsVirt)
          do k = 1, size(epsOcc)
            do j = 1, size(epsOcc)
              do i = 1, size(epsOcc)
                amplitudes(i, j, k, a, b, c) = source(i, j, k, a, b, c) / &
                                               (epsOcc(i) + epsOcc(j) + epsOcc(k) - epsVirt(a) - epsVirt(b) - epsVirt(c))
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine

end submodule
