source setBuildVars.sh -c gnu -t debug &&
cmake -B build                         &&
cmake --build build                    &&
cd build                               &&
ctest -R T_                            &&
ctest -T Coverage                      && 
lcov -c -d . -o main_coverage.info     &&
genhtml main_coverage.info --output-directory out &&
chromium out/index.html