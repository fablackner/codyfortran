 
*[main][~/projects/fProjects/codyfortranrdm]$ cd build && ctest -R T_ && cd ..
Test project /home/fabian/projects/fProjects/codyfortranrdm/build
      Start  1: T_SysInteraction_Ylm_Coulomb_BlockEq
 1/27 Test  #1: T_SysInteraction_Ylm_Coulomb_BlockEq ..........***Failed    0.40 sec
      Start  2: T_SysInteraction_Ylm_Coulomb_FullEq
 2/27 Test  #2: T_SysInteraction_Ylm_Coulomb_FullEq ...........***Failed    0.39 sec
      Start  3: T_SysInteraction_Ylm_Coulomb_StdImpl
 3/27 Test  #3: T_SysInteraction_Ylm_Coulomb_StdImpl ..........***Failed    0.39 sec
      Start  4: T_SysInteraction_Ylm_Coulomb_TwoScan
 4/27 Test  #4: T_SysInteraction_Ylm_Coulomb_TwoScan ..........***Failed    0.39 sec
      Start  5: T_SysKinetic_Ylm_Laplacian_Fedvr
 5/27 Test  #5: T_SysKinetic_Ylm_Laplacian_Fedvr ..............***Failed    0.39 sec
      Start  6: T_Utils_DerivativeFedvr
 6/27 Test  #6: T_Utils_DerivativeFedvr .......................   Passed    0.23 sec
      Start  7: T_BoseHubbard_01_ExactDiagTdciLapack
 7/27 Test  #7: T_BoseHubbard_01_ExactDiagTdciLapack ..........***Failed    0.40 sec
      Start  8: T_BoseHubbard_02_ExactDiagTdciArpack
 8/27 Test  #8: T_BoseHubbard_02_ExactDiagTdciArpack ..........***Failed    0.41 sec
      Start  9: T_BoseHubbard_03_ImagTimePropagationTdci
 9/27 Test  #9: T_BoseHubbard_03_ImagTimePropagationTdci ......***Failed    0.38 sec
      Start 10: T_BoseHubbard_04_ImagTimePropagationMctdhx
10/27 Test #10: T_BoseHubbard_04_ImagTimePropagationMctdhx ....***Failed    0.38 sec
      Start 11: T_BoseHubbard_05_ImagTimePropagationTdhx
11/27 Test #11: T_BoseHubbard_05_ImagTimePropagationTdhx ......***Failed    0.37 sec
      Start 12: T_FermiHubbard_01_ExactDiagTdciLapack
12/27 Test #12: T_FermiHubbard_01_ExactDiagTdciLapack .........***Failed    0.38 sec
      Start 13: T_FermiHubbard_02_ExactDiagTdciArpack
13/27 Test #13: T_FermiHubbard_02_ExactDiagTdciArpack .........***Failed    0.38 sec
      Start 14: T_FermiHubbard_03_ExactPropagationTdci
14/27 Test #14: T_FermiHubbard_03_ExactPropagationTdci ........***Failed    0.38 sec
      Start 15: T_FermiHubbard_04_ImagTimePropagationTdci
15/27 Test #15: T_FermiHubbard_04_ImagTimePropagationTdci .....***Failed    0.38 sec
      Start 16: T_FermiHubbard_05_ImagTimePropagationMctdhx
16/27 Test #16: T_FermiHubbard_05_ImagTimePropagationMctdhx ...***Failed    0.39 sec
      Start 17: T_FermiHubbard_06_ImagTimePropagationTdhx
17/27 Test #17: T_FermiHubbard_06_ImagTimePropagationTdhx .....***Failed    0.38 sec
      Start 18: T_H3d_01_ImagTimePropagationSingleBody
18/27 Test #18: T_H3d_01_ImagTimePropagationSingleBody ........***Failed    0.39 sec
      Start 19: T_He1d_01_ExactDiagGridBasedFullExpansion
19/27 Test #19: T_He1d_01_ExactDiagGridBasedFullExpansion .....***Failed    0.38 sec
      Start 20: T_He1d_02_ExactDiagTdci
20/27 Test #20: T_He1d_02_ExactDiagTdci .......................***Failed    0.38 sec
      Start 21: T_He1d_03_ImagTimePropagationTdhx
21/27 Test #21: T_He1d_03_ImagTimePropagationTdhx .............***Failed    0.38 sec
      Start 22: T_He1d_04_ImagTimePropagationMctdhx
22/27 Test #22: T_He1d_04_ImagTimePropagationMctdhx ...........***Failed    0.38 sec
      Start 23: T_He1d_05_TimePropagationMctdhx
23/27 Test #23: T_He1d_05_TimePropagationMctdhx ...............***Failed    0.38 sec
      Start 24: T_He1d_06_TimePropagationTd2rdm
24/27 Test #24: T_He1d_06_TimePropagationTd2rdm ...............***Failed    0.38 sec
      Start 25: T_Ne3d_01_ImagTimePropagationTdhx
25/27 Test #25: T_Ne3d_01_ImagTimePropagationTdhx .............***Failed    0.38 sec
      Start 26: T_Ne3d_02_GroundSolverTdhx
26/27 Test #26: T_Ne3d_02_GroundSolverTdhx ....................***Failed    0.39 sec
      Start 27: T_Ne3d_03_GroundSolverTdhxYlm
27/27 Test #27: T_Ne3d_03_GroundSolverTdhxYlm .................***Failed    0.40 sec

4% tests passed, 26 tests failed out of 27

Total Test time (real) =  10.27 sec

The following tests FAILED:
          1 - T_SysInteraction_Ylm_Coulomb_BlockEq (Failed)
          2 - T_SysInteraction_Ylm_Coulomb_FullEq (Failed)
          3 - T_SysInteraction_Ylm_Coulomb_StdImpl (Failed)
          4 - T_SysInteraction_Ylm_Coulomb_TwoScan (Failed)
          5 - T_SysKinetic_Ylm_Laplacian_Fedvr (Failed)
          7 - T_BoseHubbard_01_ExactDiagTdciLapack (Failed)
          8 - T_BoseHubbard_02_ExactDiagTdciArpack (Failed)
          9 - T_BoseHubbard_03_ImagTimePropagationTdci (Failed)
         10 - T_BoseHubbard_04_ImagTimePropagationMctdhx (Failed)
         11 - T_BoseHubbard_05_ImagTimePropagationTdhx (Failed)
         12 - T_FermiHubbard_01_ExactDiagTdciLapack (Failed)
         13 - T_FermiHubbard_02_ExactDiagTdciArpack (Failed)
         14 - T_FermiHubbard_03_ExactPropagationTdci (Failed)
         15 - T_FermiHubbard_04_ImagTimePropagationTdci (Failed)
         16 - T_FermiHubbard_05_ImagTimePropagationMctdhx (Failed)
         17 - T_FermiHubbard_06_ImagTimePropagationTdhx (Failed)
         18 - T_H3d_01_ImagTimePropagationSingleBody (Failed)
         19 - T_He1d_01_ExactDiagGridBasedFullExpansion (Failed)
         20 - T_He1d_02_ExactDiagTdci (Failed)
         21 - T_He1d_03_ImagTimePropagationTdhx (Failed)
         22 - T_He1d_04_ImagTimePropagationMctdhx (Failed)
         23 - T_He1d_05_TimePropagationMctdhx (Failed)
         24 - T_He1d_06_TimePropagationTd2rdm (Failed)
         25 - T_Ne3d_01_ImagTimePropagationTdhx (Failed)
         26 - T_Ne3d_02_GroundSolverTdhx (Failed)
         27 - T_Ne3d_03_GroundSolverTdhxYlm (Failed)
Errors while running CTest
Output from these tests are in: /home/fabian/projects/fProjects/codyfortranrdm/build/Testing/Temporary/LastTest.log
Use "--rerun-failed --output-on-failure" to re-run the failed cases verbosely.
*[main][~/projects/fProjects/codyfortranrdm/build]$ 


I get this result when I run the tests. i think the prob is that i have changed the folder structure in the test folder. now there are chategories in subfolders. you need to update the build and test logic to reflect the new folder structure.