! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_FieldFree) S_FieldFree

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine FieldFree_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_FieldFree_Eigen
    use M_FieldFree_Cayley

    call Say_Fabricate("fieldFree")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Nothing of an earlier fabricate survives; the variant binds what it
    ! provides
    FieldFree_Setup => NoOpProcedures_Setup
    FieldFree_ProjectBelowCut => NoOpProcedures_ApplyAbsorber
    nullify (FieldFree_Propagate, FieldFree_ApplyPhi, FieldFree_BuildMeanField, FieldFree_ApplyMeanField)
    FieldFree_cutQ = .false.
    FieldFree_meanFieldQ = .false.
    FieldFree_meanFieldBuiltQ = .false.

    !------------------------------------
    ! branch: select the variant
    !------------------------------------

    if (Json_GetExistence("fieldFree.eigen")) then
      call FieldFree_Eigen_Fabricate

    else if (Json_GetExistence("fieldFree.cayley")) then
      call FieldFree_Cayley_Fabricate

    else if (Json_GetExistence("fieldFree")) then
      error stop "fieldFree is missing one of: eigen, cayley"

    else
      ! No exact part configured → the linear integrator of the method
      ! propagates the whole one-body part
      write (*, '(1X, A, A, A)') red, 'de: FieldFree => noFieldFree', reset

    end if

  end subroutine

end submodule
