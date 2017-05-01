feelpp_add_application( poissonFE SRCS poissonFE.cpp )
feelpp_add_application( poissonCRB SRCS poissonCRB.hpp poissonCRB_main.cpp )

feelpp_add_testcase( poisson )
