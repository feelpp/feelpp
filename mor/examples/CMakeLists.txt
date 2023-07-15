OPTION(FEELPP_ENABLE_MOR_EADS "enable EADS testcase" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_HEAT1D "enable HEAT1D testcase" ON )
OPTION(FEELPP_ENABLE_MOR_HEAT2D "enable HEAT2D testcase" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_HEAT3D "enable HEAT3D testcase" ON)
OPTION(FEELPP_ENABLE_MOR_EYE2BRAIN "enable EYE2BRAIN testcase" ON)
OPTION(FEELPP_ENABLE_MOR_LINEARELASTICITY "enable LINEARELASTICITY testcase" ON)
OPTION(FEELPP_ENABLE_MOR_OPUSHEAT "enable OpusHeat testcases" ON)
OPTION(FEELPP_ENABLE_MOR_THERMALBUILDING "enable ThermalBuilding testcases" ON)
OPTION(FEELPP_ENABLE_MOR_RBHEAT "enable RBHEAT testcases" OFF)
OPTION(FEELPP_ENABLE_MOR_CABINEHEAT "enable CABINEHEAT testcases" OFF)
OPTION(FEELPP_ENABLE_MOR_FINALCABINEHEAT "enable FINALCABINEHEAT testcases" OFF)
OPTION(FEELPP_ENABLE_MOR_AD1 "enable AD1(Advection-Diffusion) testcase" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_THERMALBLOCK "enable THERMALBLOCK(many parameters) testcase" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_HELMHOLTZ "enable Helmholtz testcases" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_HEATSINK "enable Heatsink testcases" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_HEATSHIELD "enable Heatshield testcases" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_BENCHMARKGREPL "enable BenchmarkGrepl testcases" ON)
OPTION(FEELPP_ENABLE_MOR_GEOHEATDIFFUSION "enable GeoHeatDiffusion testcases" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_THERMALFIN "enable ThermalFin testcases" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_EEG "enable EEG testcase" OFF)
OPTION(FEELPP_ENABLE_MOR_CONVECTION "enable Convection testcase" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_STOKESDEIM "enable stokes with deim testcase" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_MINIMAL "enable minimal testcase" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_THERMOELECTRIC "enable thermoelectric testcase" ${FEELPP_ENABLE_MOR_ALL})
OPTION(FEELPP_ENABLE_MOR_GREPLDEIM "enable grepldeim testcase" ${FEELPP_ENABLE_MOR_ALL})


IF ( FEELPP_ENABLE_MOR_EADS )
  #add_subdirectory(eads)
ENDIF()

IF ( FEELPP_ENABLE_MOR_HEAT2D )
  add_subdirectory(heat2d)
ENDIF()

IF ( FEELPP_ENABLE_MOR_HEAT3D )
  add_subdirectory(heat3d)
ENDIF()

IF(FEELPP_ENABLE_MOR_EYE2BRAIN)
  add_subdirectory(eye2brain)
ENDIF()

IF ( FEELPP_ENABLE_MOR_HEAT1D )
  add_subdirectory(heat1d)
  add_subdirectory(unsteady_heat1d)
ENDIF()

if ( FEELPP_ENABLE_MOR_LINEARELASTICITY )
  add_subdirectory(linearelasticity3d)
endif()

if ( FEELPP_ENABLE_MOR_THERMALBUILDING )
  add_subdirectory(thermalbuilding)
endif()
if ( FEELPP_ENABLE_MOR_OPUSHEAT )
  add_subdirectory(opusheat)
endif()

IF ( FEELPP_ENABLE_MOR_RBHEAT )
  #add_subdirectory(Rbheat)
ENDIF()

IF ( FEELPP_ENABLE_MOR_CABINEHEAT )
  #add_subdirectory(Cabineheat)
ENDIF()

IF ( FEELPP_ENABLE_MOR_FINALCABINEHEAT )
  #add_subdirectory(Finalcabineheat)
ENDIF()

IF ( FEELPP_ENABLE_MOR_AD1 )
  add_subdirectory(advection-diffusion-2d-1)
ENDIF()


IF ( FEELPP_ENABLE_MOR_THERMALBLOCK )
  add_subdirectory(thermalblock)
ENDIF()


IF ( FEELPP_ENABLE_MOR_HELMHOLTZ )
add_subdirectory(helmholtz)
ENDIF()


if ( FEELPP_ENABLE_MOR_HEATSINK )
  add_subdirectory(heatsink-2d)
ENDIF()

if ( FEELPP_ENABLE_MOR_HEATSHIELD )
  add_subdirectory(heatshield)
ENDIF()

if ( FEELPP_ENABLE_MOR_BENCHMARKGREPL )
  add_subdirectory(BenchmarkGrepl)
ENDIF()

if ( FEELPP_ENABLE_MOR_GEOHEATDIFFUSION )
  add_subdirectory(GeoHeatDiffusion)
ENDIF()

if ( FEELPP_ENABLE_MOR_THERMALFIN )
  add_subdirectory(thermalfin)
ENDIF()

IF ( FEELPP_ENABLE_MOR_EEG )
  add_subdirectory(eeg)
ENDIF()

IF ( FEELPP_ENABLE_MOR_CONVECTION )
  add_subdirectory(convection)
ENDIF()

IF ( FEELPP_ENABLE_MOR_STOKESDEIM )
  add_subdirectory(stokesdeim)
ENDIF()



if ( FEELPP_ENABLE_MOR_MINIMAL )
  add_subdirectory(minimal)
endif()
if ( FEELPP_ENABLE_MOR_THERMOELECTRIC )
  add_subdirectory(thermoelectric)
endif()
if ( FEELPP_ENABLE_MOR_GREPLDEIM )
  add_subdirectory(grepldeim)
endif()
