add_custom_target(pyfeelpp)
configure_file(__init__.py.in ${CMAKE_CURRENT_BINARY_DIR}/__init__.py @ONLY)
feelpp_add_pymodule(feelpp SRCS feelpp.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(core SRCS core.cpp remotedata.cpp  DESTINATION feelpp/core )
feelpp_add_pymodule(alg SRCS alg.cpp DESTINATION feelpp/core )
if ( PETSC4PY_FOUND )
  target_compile_definitions(_alg PRIVATE -DFEELPP_HAS_PETSC4PY=1)
endif()
feelpp_add_pymodule(mesh SRCS mesh.cpp mesh1.cpp mesh2.cpp mesh3.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(discr SRCS discr.cpp discr_common.cpp discr_space_pch.cpp discr_space_pchv.cpp discr_space_pdh.cpp discr_space_pdhv.cpp discr_space_dh.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(exporter SRCS exporter.cpp hbf.cpp  DESTINATION feelpp/core )
feelpp_add_pymodule(ls SRCS ls.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(ts SRCS ts.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(forms SRCS forms.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(vf SRCS vf.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(measure SRCS measure.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(integrate SRCS integrate.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(interpolation SRCS interpolation.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(operators SRCS operators.cpp operators1.cpp operators2.cpp operators3.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(models SRCS models.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(quality SRCS quality.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(plot SRCS plot.cpp DESTINATION feelpp/core )
feelpp_add_pymodule(timing SRCS timing.cpp DESTINATION feelpp/core)
feelpp_add_pymodule(meshmover SRCS meshmover.cpp DESTINATION feelpp/core)
foreach(MODULE IN ITEMS core feelpp alg mesh discr exporter ls ts vf forms measure models integrate interpolation operators quality timing meshmover plot)
  add_dependencies(pyfeelpp _${MODULE})
endforeach()

set(PYFILES
  ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
  integrate.py
  measure.py
  quality.py
  interpolation.py
  operators.py
  timing.py
  meshmover.py
  plotly.py
  forms.py
)

feelpp_stage_python_files(
  TARGET pyfeelpp_core_pyfiles
  DESTINATION feelpp/core
  FILES ${PYFILES}
)
add_dependencies(pyfeelpp pyfeelpp_core_pyfiles)

install(FILES ${PYFILES}
        DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp/core
        COMPONENT ${FEELPP_PYTHON_INSTALL_COMPONENT})
