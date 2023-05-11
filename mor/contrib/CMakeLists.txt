include(FetchContent)


option(FEELPP_ENABLE_OPENTURNS "Enable OPENTURNS" ON)

if(FEELPP_ENABLE_OPENTURNS)
    message(STATUS "------------------------------------------------------------")
    message(STATUS "Configuring OPENTURNS")
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/openturns/cmake)

    set(USE_CERES OFF)
    set(USE_MPFR OFF)
    set(USE_MPC OFF)
    set(USE_NLOPT OFF)
    set(USE_IPOPT OFF)
    set(USE_OPENMP OFF)
    set(USE_BONMIN OFF)
    set(USE_PAGMO OFF)
    set(USE_CMINPACK OFF)
    set(BUILD_PYTHON OFF)
    set(USE_SPHINX OFF)
    set(USE_DOXYGEN OFF)
    set(BUILD_OPENTURNS_TESTS ON)
    set(OPENTURNS_HEADER_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}/")
    set(OPENTURNS_CMAKE_INSTALL_DIR "share/feelpp/openturns")
    set(OPENTURNS_DOC_INSTALL_DIR "share/doc/feelpp/openturns")
    set(FETCHCONTENT_SOURCE_DIR_OPENTURNS ${CMAKE_CURRENT_SOURCE_DIR}/openturns)

    # patch includes fix for 18.04 compilation
    FetchContent_Declare(openturns GIT_REPOSITORY https://github.com/feelpp/openturns.git GIT_TAG 4198b836e5fecc0757241f2dc74ec5cc6f79fe40)

    FetchContent_MakeAvailable(openturns)
    add_library(feelpp_openturns INTERFACE)
    set_target_properties(OT PROPERTIES OUTPUT_NAME "feelpp_openturns")
    target_link_libraries(feelpp_openturns INTERFACE OT)
# #    target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_OPENTURNS)
    add_library(Feelpp::feelpp_openturns ALIAS feelpp_openturns)
    target_include_directories(feelpp_openturns INTERFACE
        $<BUILD_INTERFACE:${openturns_SOURCE_DIR}/lib/src>
        $<BUILD_INTERFACE:${openturns_SOURCE_DIR}/lib/src/Base>
        $<BUILD_INTERFACE:${openturns_SOURCE_DIR}/lib/include/>
        $<BUILD_INTERFACE:${openturns_BINARY_DIR}/lib/include/>
        $<INSTALL_INTERFACE:include/feelpp>)

    # Use the file(GLOB ...) command to list directories
    file(GLOB BASE_DIRS RELATIVE ${openturns_SOURCE_DIR}/lib/src/Base ${openturns_SOURCE_DIR}/lib/src/Base/*)

    foreach(dir in 
        python/src/
        python/doc/themes/
        lib/src/
        lib/src/Base/Common/
        lib/src/Base/
        lib/src/Base/Diff/
        lib/src/Base/MetaModel/
        lib/src/Base/Algo/
        lib/src/Base/Graph/
        lib/src/Base/Geom/
        lib/src/Base/Experiments/
        lib/src/Base/Solver/
        lib/src/Base/Func/
        lib/src/Base/Func/SpecFunc/
        lib/src/Base/Type/
        lib/src/Base/IterativeStat/
        lib/src/Base/Stat/
        lib/src/Base/Optim/
        lib/src/Uncertainty/
        lib/src/Uncertainty/Algorithm/Analytical/
        lib/src/Uncertainty/Algorithm/Transformation/
        lib/src/Uncertainty/Algorithm/Transformation/ProcessTransformation/
        lib/src/Uncertainty/Algorithm/Transformation/MarginalTransformation/
        lib/src/Uncertainty/Algorithm/Transformation/NatafEllipticalDistribution/
        lib/src/Uncertainty/Algorithm/Transformation/Rosenblatt/
        lib/src/Uncertainty/Algorithm/Transformation/NatafEllipticalCopula/
        lib/src/Uncertainty/Algorithm/Transformation/NatafIndependentCopula/
        lib/src/Uncertainty/Algorithm/OrthogonalBasis/
        lib/src/Uncertainty/Algorithm/
        lib/src/Uncertainty/Algorithm/Classification/
        lib/src/Uncertainty/Algorithm/MetaModel/
        lib/src/Uncertainty/Algorithm/MetaModel/LinearModel/
        lib/src/Uncertainty/Algorithm/MetaModel/Kriging/
        lib/src/Uncertainty/Algorithm/MetaModel/FunctionalChaos/
        lib/src/Uncertainty/Algorithm/Simulation/
        lib/src/Uncertainty/Algorithm/WeightedExperiments/
        lib/src/Uncertainty/Algorithm/Sensitivity/
        lib/src/Uncertainty/Algorithm/Optimization/
        lib/src/Uncertainty/Distribution/
        lib/src/Uncertainty/StatTests/
        lib/src/Uncertainty/Process/
        lib/src/Uncertainty/Model/
        lib/src/Uncertainty/Bayesian/
        lib/include/)
            target_include_directories(feelpp_openturns INTERFACE
                $<BUILD_INTERFACE:${openturns_SOURCE_DIR}/${dir}>)
        endforeach()

        install(TARGETS feelpp_openturns OT EXPORT feelpp-mor-export-targets
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )
        endif()