

set(PYFILES test_core.py test_remotedata.py test_discr.py test_exporter.py test_mesh.py )
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)