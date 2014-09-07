## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "Feel++")
set(CTEST_NIGHTLY_START_TIME "00:00:00 CEST")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=Feel%2B%2B")
set(CTEST_DROP_SITE_CDASH TRUE)

set(CTEST_UPDATE_TYPE "git")

set(CTEST_PROJECT_SUBPROJECTS
feelpp
testcore
testfilters
testmaterial
testleaks
testcrb
testdiscr
testmesh
testalg
testinterpolation
testintegration
testpoly
testts
testvf
testsuite
benchmarks
doc
crb
)
