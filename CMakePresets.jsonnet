// CMakePresets.jsonnet - Template for generating CMakePresets.json
// Generate with: ./scripts/generate-presets.sh
// Or manually: jsonnet CMakePresets.jsonnet | python3 -m json.tool --indent 3 --no-ensure-ascii > CMakePresets.json

local compilers = ['gcc', 'clang'];
local cppStds = ['20', '23'];
local buildTypes = {
  release: 'Release',
  debug: 'Debug',
  relwithdebinfo: 'RelWithDebInfo',
};
local components = ['feelpp', 'testsuite', 'quickstart', 'toolboxes', 'mor', 'python'];
local packageManagers = ['cmake', 'spack', 'conan', 'vcpkg'];
local spackEnvironments = ['cpu/openmpi', 'cpu/openmpi5'];
local spackPresetName(env) = 'release-clang-spack-' + std.strReplace(env, '/', '-');

// ============================================================================
// Helper Functions
// ============================================================================

local capitalize(str) = 
  std.asciiUpper(str[0:1]) + str[1:];

local componentDisplayName(comp) = {
  feelpp: 'Feel++ Library',
  testsuite: 'Test Suite',
  quickstart: 'Quickstart',
  toolboxes: 'Toolboxes',
  mor: 'MOR (Model Order Reduction)',
  python: 'Python Bindings',
}[comp];

// ============================================================================
// Base Presets (Foundation)
// ============================================================================

local defaultPreset = {
  name: 'default',
  displayName: 'Default Config',
  description: 'Default build using Ninja generator',
  generator: 'Ninja',
  binaryDir: '${sourceDir}/build/${presetName}$env{DISTRIBUTION}$env{ARCH}',
  cacheVariables: {
    CMAKE_INSTALL_PREFIX: '${sourceDir}/install/${presetName}$env{DISTRIBUTION}$env{ARCH}',
    CMAKE_CXX_COMPILER: 'clang++',
    CMAKE_C_COMPILER: 'clang',
    FEELPP_STD_CPP: '23',  // Default to C++23
    FEELPP_ENABLE_TOOLBOXES: 'ON',
    FEELPP_ENABLE_MOR: 'ON',
    FEELPP_ENABLE_FEELPP_PYTHON: 'ON',
    FEELPP_ENABLE_FMILIB: 'OFF',
    FEELPP_ENABLE_BENCHMARKS: 'ON',
    FEELPP_USE_EXTERNAL_EIGEN3: 'OFF',
    FEELPP_USE_EXTERNAL_PYBIND11: 'ON',
    CMAKE_EXPORT_COMPILE_COMMANDS: 'TRUE',
    CMAKE_POLICY_VERSION_MINIMUM: '3.5',
    // Enable the maintained toolbox set in the main build.
    FEELPP_TOOLBOXES_ENABLE_HEAT: 'ON',
    FEELPP_TOOLBOXES_ENABLE_CFPDE: 'ON',
    FEELPP_TOOLBOXES_ENABLE_ELECTRIC: 'ON',
    FEELPP_TOOLBOXES_ENABLE_THERMOELECTRIC: 'ON',
    FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS: 'ON',
    FEELPP_TOOLBOXES_ENABLE_SOLIDMECHANICS: 'ON',
    FEELPP_TOOLBOXES_ENABLE_FSI: 'OFF',
    FEELPP_TOOLBOXES_ENABLE_ADVECTION: 'OFF',
    FEELPP_TOOLBOXES_ENABLE_LEVELSET: 'OFF',
    FEELPP_TOOLBOXES_ENABLE_MULTIFLUID: 'OFF',
    FEELPP_TOOLBOXES_ENABLE_HDG: 'OFF',
    FEELPP_TOOLBOXES_ENABLE_MAXWELL: 'OFF',
  },
  environment: {
    LDFLAGS: '-Wl,--copy-dt-needed-entries',
    NINJA_STATUS: '[%f/%t %p%%] %r⚙ %o/s %es | ',
  },
  vendor: {
    'example.com/ExampleIDE/1.0': {
      autoFormat: true,
    },
  },
};

// ============================================================================
// Hidden Mixin Presets
// ============================================================================

local warningsPreset = {
  name: 'warnings',
  hidden: true,
  description: 'Enable compiler warnings on linux and osx',
  cacheVariables: {
    CMAKE_CXX_FLAGS: '-Wall -Wextra -Wpedantic -Wconversion -Wunused',
    CMAKE_C_FLAGS: '-Wall -Wextra -Wpedantic -Wconversion -Wunused',
  },
};

local usrlocalPreset = {
  name: 'usrlocal',
  hidden: true,
  description: 'Install in /usr/local',
  cacheVariables: {
    CMAKE_INSTALL_PREFIX: '/usr/local',
  },
};

local macosxPreset = {
  name: 'macosx',
  hidden: true,
  description: 'macOS build settings',
  environment: {
    LDFLAGS: '',
  },
  cacheVariables: {
    FEELPP_ENABLE_FFTW: 'ON',
    FEELPP_ENABLE_GLPK: 'ON',
  },
};

// C++ Standard presets
local cppStdPreset(std) = {
  name: 'cpp' + std,
  hidden: true,
  description: 'Enable C++' + std,
  cacheVariables: {
    FEELPP_STD_CPP: std,
  },
};

// Compiler presets
local compilerPreset(compiler) = {
  name: compiler,
  hidden: true,
  description: 'Use ' + compiler + ' compiler',
  cacheVariables: {
    CMAKE_CXX_COMPILER: if compiler == 'gcc' then 'g++' else compiler + '++',
    CMAKE_C_COMPILER: compiler,
  },
};

// Package manager presets
local spackPreset = {
  name: 'spack',
  hidden: true,
  displayName: 'spack package manager',
  description: 'spack config',
  inherits: ['default'],
  cacheVariables: {
    CMAKE_INSTALL_RPATH_USE_LINK_PATH: 'ON',
    FEELPP_ALLOW_AMBIENT_PYTHON: 'ON',
    FEELPP_USE_EXTERNAL_CLN: 'ON',
    FEELPP_USE_EXTERNAL_EIGEN3: 'ON',
    FEELPP_ENABLE_VTK: 'OFF',
    USE_VTK: 'OFF',
    FEELPP_ENABLE_OPENTURNS: 'OFF',
    PYBIND11_FINDPYTHON: 'ON',
    CPR_USE_SYSTEM_CURL: 'ON',
    CPR_USE_SYSTEM_LIB_PSL: 'OFF',
    CURL_NO_CURL_CMAKE: 'ON',
    FEELPP_ENABLE_OMC: 'OFF',
  },
  environment: {
    VERBOSE: '1',
  },
};

// ============================================================================
// Build Type Base Presets
// ============================================================================

local releasePreset = {
  name: 'release',
  hidden: false,
  displayName: 'Release | no package manager',
  description: 'Plain release build without dependency handling.',
  inherits: ['default'],
  cacheVariables: {
    CMAKE_BUILD_TYPE: 'Release',
  },
};

local debugPreset = {
  name: 'debug',
  hidden: false,
  displayName: 'Debug | no package manager',
  description: 'Plain debug build without dependency handling.',
  inherits: ['default'],
  cacheVariables: {
    CMAKE_BUILD_TYPE: 'Debug',
    CMAKE_VERBOSE_MAKEFILE: 'TRUE',
    CMAKE_MESSAGE_LOG_LEVEL: 'VERBOSE',
    H5PP_ENABLE_TESTS: 'TRUE',
    H5PP_BUILD_EXAMPLES: 'TRUE',
    H5PP_ENABLE_ASAN: 'TRUE',
    H5PP_ENABLE_PCH: 'FALSE',
    H5PP_ENABLE_CCACHE: 'FALSE',
    CMAKE_INTERPROCEDURAL_OPTIMIZATION: 'FALSE',
    CMAKE_COMPILE_WARNING_AS_ERROR: 'FALSE',
    FEELPP_ENABLE_TOOLBOXES: 'ON',
    FEELPP_ENABLE_MOR: 'OFF',
    FEELPP_ENABLE_FEELPP_PYTHON: 'OFF',
    FEELPP_INSTANTIATION_ORDER_MAX: '1',
    FEELPP_MESH_MAX_ORDER: '1',
  },
  environment: {
    VERBOSE: '1',
  },
};

local devPreset = {
  name: 'dev',
  hidden: false,
  displayName: 'Development | no package manager',
  description: 'development config',
  inherits: ['default'],
  cacheVariables: {
    CMAKE_BUILD_TYPE: 'Release',
    CMAKE_VERBOSE_MAKEFILE: 'TRUE',
    CMAKE_MESSAGE_LOG_LEVEL: 'VERBOSE',
    H5PP_ENABLE_TESTS: 'TRUE',
    H5PP_BUILD_EXAMPLES: 'TRUE',
    H5PP_ENABLE_ASAN: 'TRUE',
    H5PP_ENABLE_PCH: 'FALSE',
    H5PP_ENABLE_CCACHE: 'FALSE',
    CMAKE_INTERPROCEDURAL_OPTIMIZATION: 'FALSE',
    CMAKE_COMPILE_WARNING_AS_ERROR: 'FALSE',
    FEELPP_ENABLE_TOOLBOXES: 'ON',
    FEELPP_ENABLE_MOR: 'ON',
    FEELPP_ENABLE_FEELPP_PYTHON: 'ON',
    FEELPP_INSTANTIATION_ORDER_MAX: '1',
    FEELPP_MESH_MAX_ORDER: '1',
  },
  environment: {
    VERBOSE: '1',
  },
};

// C++23 default preset
local defaultCpp23Preset = {
  name: 'default-cpp23',
  hidden: false,
  displayName: 'C++23 | no package manager',
  description: 'Default build with C++23 enabled',
  inherits: ['cpp23', 'default'],
  cacheVariables: {
    CMAKE_BUILD_TYPE: 'Release',
  },
};

// Special debug preset (dbg)
local dbgPreset = {
  name: 'dbg',
  displayName: 'Debug Config',
  description: 'Debug build using Ninja generator',
  generator: 'Ninja',
  binaryDir: '${sourceDir}/build/dbg',
  cacheVariables: {
    CMAKE_CXX_COMPILER: 'clang++',
    CMAKE_BUILD_TYPE: 'Debug',
    FEELPP_ENABLE_TOOLBOXES: {
      type: 'BOOL',
      value: 'ON',
    },
    FEELPP_ENABLE_MOR: 'OFF',
    FEELPP_INSTANTIATION_ORDER_MAX: '1',
    FEELPP_MESH_MAX_ORDER: '1',
  },
  environment: {
    MY_ENVIRONMENT_VARIABLE: 'Test',
  },
  vendor: {
    'example.com/ExampleIDE/1.0': {
      autoFormat: true,
    },
  },
};

// ============================================================================
// Matrix Generation: Build Type + Package Manager + Compiler + C++ Std
// ============================================================================

// Simple build-type + cmake presets
local buildTypeCmakePreset(buildType) = {
  name: buildType + '-cmake',
  displayName: capitalize(buildType) + ' | cmake package manager',
  description: 'Uses a custom wrapper for external_project_add at CMake configure time',
  inherits: [buildType, 'default'],
};

// Build-type + cpp-std + cmake
local buildTypeCppStdCmakePreset(buildType, cppStd) = {
  name: buildType + '-cpp' + cppStd + '-cmake',
  displayName: capitalize(buildType) + ' | C++' + cppStd + ' | cmake package manager',
  description: 'Uses a custom wrapper for external_project_add at CMake configure time',
  inherits: [buildType, 'cpp' + cppStd, 'default'],
};

// Build-type + compiler + cmake
local buildTypeCompilerCmakePreset(buildType, compiler) = {
  name: buildType + '-' + compiler + '-cmake',
  displayName: capitalize(buildType) + ' | ' + compiler + ' | cmake package manager',
  inherits: [compiler, buildType + '-cmake'],
};

// Build-type + compiler + dev + cmake (special for clang dev)
local buildTypeCompilerDevCmakePreset(buildType, compiler) = {
  name: buildType + '-' + compiler + '-cmake-dev',
  displayName: capitalize(buildType) + ' | ' + compiler + ' | cmake package manager',
  inherits: ['dev', compiler, buildType + '-cmake'],
};

// Build-type + compiler + cpp-std + cmake
local buildTypeCompilerCppStdCmakePreset(buildType, compiler, cppStd) = {
  name: buildType + '-' + compiler + '-cpp' + cppStd + '-cmake',
  displayName: capitalize(buildType) + ' | ' + compiler + ' | C++' + cppStd + ' | cmake package manager',
  inherits: [compiler, 'cpp' + cppStd, buildType + '-cmake'],
};

// Build-type + compiler + spack
local buildTypeCompilerSpackPreset(buildType, compiler) = {
  name: buildType + '-' + compiler + '-spack',
  displayName: capitalize(buildType) + ' | ' + compiler + ' | spack package manager',
  inherits: ['spack', compiler, buildType + '-cmake'],
};

// Build-type + compiler + spack + macOS
local buildTypeCompilerSpackMacosxPreset(buildType, compiler) = {
  name: buildType + '-' + compiler + '-spack-macosx',
  displayName: capitalize(buildType) + ' | ' + compiler + ' | spack package manager | macOS',
  inherits: ['macosx', 'spack', compiler, buildType + '-cmake'],
  environment: {
    LDFLAGS: '',
  },
};

// Build-type + compiler + cpp-std + spack
local buildTypeCompilerCppStdSpackPreset(buildType, compiler, cppStd) = {
  name: buildType + '-' + compiler + '-cpp' + cppStd + '-spack',
  displayName: capitalize(buildType) + ' | ' + compiler + ' | C++' + cppStd + ' | spack package manager',
  inherits: ['cpp' + cppStd, 'spack', compiler, buildType + '-cmake'],
};

// ============================================================================
// Component Presets
// ============================================================================

local componentCacheVars = {
  feelpp: {
    FEELPP_COMPONENT: 'feelpp',
    FEELPP_ENABLE_MOR: 'OFF',
    FEELPP_ENABLE_TOOLBOXES: 'OFF',
    FEELPP_ENABLE_FEELPP_PYTHON: 'OFF',
    FEELPP_ENABLE_PYTHON: 'ON',
    FEELPP_ENABLE_TESTS: 'OFF',
    FEELPP_ENABLE_FMILIB: 'OFF',
    FEELPP_ENABLE_BENCHMARKS: 'OFF',
    FEELPP_ENABLE_QUICKSTART: 'OFF',
    FEELPP_USE_EXTERNAL_CLN: 'ON',
  },
  toolboxes: {
    FEELPP_COMPONENT: 'toolboxes',
    FEELPP_ENABLE_FEELPP_PYTHON: 'OFF',
    FEELPP_ENABLE_PYTHON: 'ON',
    FEELPP_TOOLBOXES_ENABLE_PYTHON: 'ON',
    // Enable the maintained toolbox set for the toolbox component build.
    FEELPP_TOOLBOXES_ENABLE_HEAT: 'ON',
    FEELPP_TOOLBOXES_ENABLE_CFPDE: 'ON',
    FEELPP_TOOLBOXES_ENABLE_ELECTRIC: 'ON',
    FEELPP_TOOLBOXES_ENABLE_THERMOELECTRIC: 'ON',
    FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS: 'ON',
    FEELPP_TOOLBOXES_ENABLE_SOLIDMECHANICS: 'ON',
    FEELPP_TOOLBOXES_ENABLE_FSI: 'OFF',
    FEELPP_TOOLBOXES_ENABLE_ADVECTION: 'OFF',
    FEELPP_TOOLBOXES_ENABLE_LEVELSET: 'OFF',
    FEELPP_TOOLBOXES_ENABLE_MULTIFLUID: 'OFF',
    FEELPP_TOOLBOXES_ENABLE_HDG: 'OFF',
    FEELPP_TOOLBOXES_ENABLE_MAXWELL: 'OFF',
  },
  mor: {
    FEELPP_COMPONENT: 'mor',
    FEELPP_ENABLE_RESEARCH: 'OFF',
    FEELPP_ENABLE_OPENTURNS: 'ON',
    FEELPP_ENABLE_FEELPP_PYTHON: 'OFF',
    FEELPP_ENABLE_PYTHON: 'ON',
    FEELPP_MOR_ENABLE_PYTHON: 'ON',
  },
  python: {
    FEELPP_COMPONENT: 'python',
  },
  testsuite: {
    FEELPP_COMPONENT: 'testsuite',
  },
  quickstart: {
    FEELPP_COMPONENT: 'quickstart',
  },
};

local componentPreset(component) = {
  name: component,
  inherits: ['clang', 'release-cmake'],
  displayName: component + ' | clang | release | cmake package manager',
  description: 'Build only the ' + componentDisplayName(component) + ' component',
  cacheVariables: componentCacheVars[component],
};

// Special component presets
local feelppUsrlocalPreset = {
  name: 'feelpp-usrlocal',
  inherits: ['usrlocal', 'feelpp'],
  displayName: 'feelpp | usrlocal | clang | release | cmake package manager',
  description: 'Build only the Feel++ library Component and install in /usr/local',
};

local feelppCpp20SpackPreset = {
  name: 'feelpp-cpp20-spack',
  inherits: ['feelpp', 'cpp20', 'spack'],
  displayName: 'feelpp | clang | cpp20 | release | spack package manager',
  description: 'Build only the Feel++ library Component',
};

local feelppCpp23SpackPreset = {
  name: 'feelpp-cpp23-spack',
  inherits: ['feelpp', 'cpp23', 'spack'],
  displayName: 'feelpp | clang | cpp23 | release | spack package manager',
  description: 'Build only the Feel++ library Component',
};

local feelppSpecxPreset = {
  name: 'feelpp+specx',
  inherits: ['clang', 'release-cmake'],
  displayName: 'feelpp+specx | clang | release | cmake package manager',
  description: 'Build only the Feel++ library Component with SpecX',
  cacheVariables: {
    FEELPP_COMPONENT: 'feelpp',
    FEELPP_ENABLE_MOR: 'OFF',
    FEELPP_ENABLE_TOOLBOXES: 'OFF',
    FEELPP_ENABLE_FEELPP_PYTHON: 'OFF',
    FEELPP_ENABLE_TESTS: 'OFF',
    FEELPP_ENABLE_FMILIB: 'OFF',
    FEELPP_ENABLE_BENCHMARKS: 'ON',
    FEELPP_ENABLE_SPECX: 'ON',
  },
};

// Special debug component presets
local morDbgPreset = {
  name: 'mor-dbg',
  inherits: 'mor',
  displayName: 'mor-dbg',
  description: 'Build only the Feel++ MOR Component in Debug',
  binaryDir: '${sourceDir}/build/mor-dbg',
  cacheVariables: {
    CMAKE_CXX_FLAGS_DEBUG: '-g -O0',
    CMAKE_INSTALL_PREFIX: '${sourceDir}/build/mor-dbg/install/',
    FEELPP_COMPONENT: 'mor',
    FEELPP_ENABLE_RESEARCH: 'OFF',
    CMAKE_BUILD_TYPE: 'Debug',
    FEELPP_ENABLE_OPENTURNS: 'ON',
    CMAKE_VERBOSE_MAKEFILE: 'OFF',
  },
};

local researchMorPreset = {
  name: 'research-mor',
  inherits: 'default',
  displayName: 'research-mor',
  description: 'Build only the Feel++ MOR and Research Components',
  binaryDir: '${sourceDir}/build/research-mor',
  cacheVariables: {
    FEELPP_COMPONENT: 'mor',
    FEELPP_ENABLE_RESEARCH: 'ON',
  },
};

local feelppPythonPreset = {
  name: 'feelpp-python',
  inherits: ['clang', 'release-cmake'],
  displayName: 'python | clang | release | cmake package manager',
  description: 'Build only the Feel++ python Component',
  cacheVariables: {
    FEELPP_COMPONENT: 'python',
  },
};

local feelppPythonDbgPreset = {
  name: 'feelpp-python-dbg',
  inherits: 'default',
  displayName: 'feelpp-python-dbg',
  description: 'Build only the Feel++ python Component in Debug mode',
  binaryDir: '${sourceDir}/build/python-dbg',
  cacheVariables: {
    CMAKE_CXX_FLAGS_DEBUG: '-g -O0',
    CMAKE_INSTALL_PREFIX: '${sourceDir}/build/python-dbg/install/',
    FEELPP_COMPONENT: 'python',
    CMAKE_BUILD_TYPE: 'Debug',
    CMAKE_VERBOSE_MAKEFILE: 'OFF',
  },
};

local morPythonPreset = {
  name: 'mor_python',
  inherits: 'default',
  displayName: 'mor_python',
  description: 'Build only the Feel++ MOR and Python Components',
  binaryDir: '${sourceDir}/build/mor_python',
  cacheVariables: {
    FEELPP_COMPONENT: 'mor;python',
  },
};

local doxPreset = {
  name: 'doxygen',
  inherits: 'default',
  displayName: 'doxygen',
  description: 'Doxygen config',
  binaryDir: '${sourceDir}/build/doxygen',
  cacheVariables: {
    FEELPP_ENABLE_DOXYGEN: 'ON',
  },
};

local focusedToolboxCacheVars = {
  FEELPP_TOOLBOXES_ENABLE_ADVECTION: 'OFF',
  FEELPP_TOOLBOXES_ENABLE_CFPDE: 'OFF',
  FEELPP_TOOLBOXES_ENABLE_ELECTRIC: 'OFF',
  FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS: 'OFF',
  FEELPP_TOOLBOXES_ENABLE_FSI: 'OFF',
  FEELPP_TOOLBOXES_ENABLE_HDG: 'OFF',
  FEELPP_TOOLBOXES_ENABLE_HEAT: 'OFF',
  FEELPP_TOOLBOXES_ENABLE_LEVELSET: 'OFF',
  FEELPP_TOOLBOXES_ENABLE_MAXWELL: 'OFF',
  FEELPP_TOOLBOXES_ENABLE_MULTIFLUID: 'OFF',
  FEELPP_TOOLBOXES_ENABLE_SOLIDMECHANICS: 'OFF',
  FEELPP_TOOLBOXES_ENABLE_THERMOELECTRIC: 'OFF',
};

local pythonCiBaseCacheVars = {
  FEELPP_ENABLE_ASCIIDOCTOR: 'OFF',
  FEELPP_ENABLE_BENCHMARKS: 'OFF',
  FEELPP_ENABLE_DOCUMENTATION: 'OFF',
  FEELPP_ENABLE_QUICKSTART: 'OFF',
  FEELPP_ENABLE_RESEARCH: 'OFF',
  FEELPP_ENABLE_TESTS: 'OFF',
};

local pythonCiPreset(name, description, displayName, extraCacheVariables={}) = {
  name: name,
  inherits: 'default',
  description: description,
  displayName: displayName,
  binaryDir: '${sourceDir}/build/' + name,
  cacheVariables: pythonCiBaseCacheVars + {
    CMAKE_INSTALL_PREFIX: '${sourceDir}/install/' + name,
  } + extraCacheVariables,
};

local pythonCiConfigurePresets = [
  pythonCiPreset(
    'python-core-ci',
    'Fresh build directory for core Python CI validation',
    'python core ci',
    {
      FEELPP_ENABLE_MOR: 'OFF',
      FEELPP_ENABLE_TOOLBOXES: 'OFF',
    }
  ),
  pythonCiPreset(
    'python-electric-ci',
    'Fresh build directory for electric toolbox Python CI validation',
    'python electric ci',
    {
      FEELPP_ENABLE_MOR: 'OFF',
    } + focusedToolboxCacheVars + {
      FEELPP_TOOLBOXES_ENABLE_ELECTRIC: 'ON',
    }
  ),
  pythonCiPreset(
    'python-mor-ci',
    'Fresh build directory for MOR Python CI validation',
    'python mor ci',
    focusedToolboxCacheVars + {
      FEELPP_TOOLBOXES_ENABLE_ELECTRIC: 'ON',
    }
  ),
  pythonCiPreset(
    'python-heat-ci',
    'Fresh build directory for heat toolbox Python CI validation',
    'python heat ci',
    {
      FEELPP_ENABLE_MOR: 'OFF',
    } + focusedToolboxCacheVars + {
      FEELPP_TOOLBOXES_ENABLE_HEAT: 'ON',
    }
  ),
  pythonCiPreset(
    'python-solid-ci',
    'Fresh build directory for solid toolbox Python CI validation',
    'python solid ci',
    {
      FEELPP_ENABLE_MOR: 'OFF',
    } + focusedToolboxCacheVars + {
      FEELPP_TOOLBOXES_ENABLE_SOLIDMECHANICS: 'ON',
    }
  ),
  pythonCiPreset(
    'python-hdg-ci',
    'Fresh build directory for HDG toolbox Python CI validation',
    'python hdg ci',
    {
      FEELPP_ENABLE_MOR: 'OFF',
    } + focusedToolboxCacheVars + {
      FEELPP_TOOLBOXES_ENABLE_HDG: 'ON',
    }
  ),
];

local pythonCiBuildPresetNames = [
  'python-core-ci',
  'python-electric-ci',
  'python-mor-ci',
  'python-heat-ci',
  'python-solid-ci',
  'python-hdg-ci',
];

// ============================================================================
// Aggregate All Configure Presets
// ============================================================================

local configurePresets = 
  // Base and mixins
  [
    defaultPreset,
    warningsPreset,
    usrlocalPreset,
    macosxPreset,
  ] +
  // C++ standard presets
  [cppStdPreset(std) for std in cppStds] +
  // Compiler presets
  [compilerPreset(c) for c in compilers] +
  // Build types
  [
    releasePreset,
    debugPreset,
    devPreset,
    defaultCpp23Preset,
    dbgPreset,
  ] +
  // Package manager presets
  [spackPreset] +
  // Build type + cmake
  [buildTypeCmakePreset('release'), buildTypeCmakePreset('debug')] +
  // Build type + cpp-std + cmake
  [buildTypeCppStdCmakePreset('release', std) for std in cppStds] +
  // Build type + compiler + cmake (with dev variant for clang)
  [buildTypeCompilerCmakePreset('release', 'gcc')] +
  [buildTypeCompilerDevCmakePreset('release', 'clang')] +
  [buildTypeCompilerCmakePreset('release', 'clang')] +
  [buildTypeCompilerCmakePreset('debug', c) for c in compilers] +
  // Build type + compiler + cpp-std + cmake
  std.flattenArrays([
    [buildTypeCompilerCppStdCmakePreset('release', c, std) for std in cppStds]
    for c in compilers
  ]) +
  // Build type + compiler + spack
  [buildTypeCompilerSpackPreset('release', 'clang')] +
  [{ name: spackPresetName(env), inherits: ['release-clang-spack'] } for env in spackEnvironments] +
  [buildTypeCompilerSpackMacosxPreset('release', 'clang')] +
  [buildTypeCompilerSpackPreset('debug', 'clang')] +
  [buildTypeCompilerSpackMacosxPreset('debug', 'clang')] +
  // Build type + compiler + cpp-std + spack
  [buildTypeCompilerCppStdSpackPreset('release', 'clang', '20')] +
  // Component presets
  [componentPreset(comp) for comp in components] +
  // Special component presets
  [
    feelppUsrlocalPreset,
    feelppCpp20SpackPreset,
    feelppCpp23SpackPreset,
    feelppSpecxPreset,
    morDbgPreset,
    researchMorPreset,
    feelppPythonPreset,
    feelppPythonDbgPreset,
    morPythonPreset,
    doxPreset,
  ] +
  pythonCiConfigurePresets;

// ============================================================================
// Build Presets
// ============================================================================

local buildPreset(configName, jobs=25) = {
  name: configName,
  configurePreset: configName,
  jobs: jobs,
};

local buildPresets = [
  buildPreset('default'),
  buildPreset('release'),
  buildPreset('release-cmake'),
  buildPreset('debug'),
  buildPreset('debug-cmake'),
  buildPreset('dbg', 2),
  buildPreset('doxygen'),
] +
// Build type + cpp-std + cmake
[buildPreset('release-cpp' + std + '-cmake') for std in cppStds] +
// Build type + compiler + cmake
[buildPreset('release-' + c + '-cmake') for c in compilers] +
[buildPreset('debug-' + c + '-cmake') for c in compilers] +
// Build type + compiler + cpp-std + cmake (matrix)
std.flattenArrays([
  [buildPreset('release-' + c + '-cpp' + std + '-cmake') for std in cppStds]
  for c in compilers
]) +
// Build type + compiler + spack
[buildPreset('release-clang-spack')] +
[buildPreset(spackPresetName(env)) for env in spackEnvironments] +
[buildPreset('release-clang-spack-macosx')] +
[buildPreset('debug-clang-spack')] +
[buildPreset('debug-clang-spack-macosx')] +
[buildPreset('release-clang-cpp20-spack')] +
// Component presets
[buildPreset(comp) for comp in components] +
// Special presets
[
  buildPreset('feelpp-usrlocal'),
  buildPreset('feelpp-cpp20-spack'),
  buildPreset('feelpp-cpp23-spack'),
  buildPreset('feelpp+specx'),
  buildPreset('mor-dbg'),
  buildPreset('research-mor'),
  buildPreset('feelpp-python'),
  buildPreset('feelpp-python-dbg'),
  buildPreset('mor_python'),
] +
[buildPreset(name, 20) for name in pythonCiBuildPresetNames];

// ============================================================================
// Test Presets
// ============================================================================

local testPreset(configName, extraConfig={}) = {
  name: configName,
  configurePreset: configName,
  output: { outputOnFailure: true },
} + extraConfig;

// Test preset with retry for failed tests (rerun up to 3 times before giving up)
local testPresetWithRetry(configName, extraConfig={}) = {
  name: configName,
  configurePreset: configName,
  output: { outputOnFailure: true },
  execution: {
    repeat: {
      mode: 'until-pass',
      count: 3,
    },
  },
} + extraConfig;

local pythonCiTestPresets = [
  testPreset('python-core-ci', {
    execution: { jobs: 4 },
  }),
  testPreset('python-electric-ci', {
    execution: { jobs: 4 },
    filter: {
      include: {
        name: '^(feelpp_qs_python-(mpi-)?core|feelpp_toolbox_python-tests-(mpi-)?electric)$',
      },
    },
  }),
  testPreset('python-mor-ci', {
    execution: { jobs: 4 },
    filter: {
      include: {
        name: '^(feelpp_qs_python-(mpi-)?core|feelpp_mor_python-tests-(mpi-)?mor)$',
      },
    },
  }),
  testPreset('python-heat-ci', {
    execution: { jobs: 4 },
    filter: {
      include: {
        name: '^(feelpp_qs_python-(mpi-)?core|feelpp_toolbox_python-tests-(mpi-)?heat(-thermo2d|-time-stepping)?|feelpp_toolbox_python-tests-(mpi-)?interpolation)$',
      },
    },
  }),
  testPreset('python-solid-ci', {
    execution: { jobs: 4 },
    filter: {
      include: {
        name: '^(feelpp_qs_python-(mpi-)?core|feelpp_toolbox_python-tests-(mpi-)?solid)$',
      },
    },
  }),
  testPreset('python-hdg-ci', {
    execution: { jobs: 4 },
    filter: {
      include: {
        name: '^(feelpp_qs_python-(mpi-)?core|feelpp_toolbox_python-tests-(mpi-)?hdg)$',
      },
    },
  }),
];

// ============================================================================
// Workflow Presets (CMake 3.25+)
// ============================================================================
// Workflow presets define a sequence of steps (configure, build, test)
// that can be run with: cmake --workflow --preset <name>

local workflowPreset(name, configPreset, buildPreset, testPreset) = {
  name: name,
  displayName: 'Workflow: ' + name,
  description: 'Configure, build, and test ' + name,
  steps: [
    { type: 'configure', name: configPreset },
    { type: 'build', name: buildPreset },
    { type: 'test', name: testPreset },
  ],
};

// Simple workflow for presets where all names match
local simpleWorkflow(name) = workflowPreset(name, name, name, name);

// Docker workflow - used by Dockerfiles for build/test cycle
local dockerWorkflow(name) = {
  name: name + '-docker',
  displayName: 'Docker: ' + name,
  description: 'Configure, build, and test ' + name + ' (for Docker builds)',
  steps: [
    { type: 'configure', name: name },
    { type: 'build', name: name },
    { type: 'test', name: name },
  ],
};

local workflowPresets =
  // Component workflows (feelpp, testsuite, quickstart, toolboxes, mor, python)
  [simpleWorkflow(comp) for comp in components] +
  // Docker workflows for CI
  [dockerWorkflow(comp) for comp in components] +
  // Default/full build workflow
  [simpleWorkflow('default')] +
  [dockerWorkflow('default')] +
  // Release workflows
  [simpleWorkflow('release')] +
  [simpleWorkflow('release-cmake')] +
  // Debug workflows
  [simpleWorkflow('debug')] +
  [simpleWorkflow('debug-cmake')];

local testPresets = [
  testPreset('default', { execution: { jobs: 4 } }),
  testPreset('release', { inherits: 'default' }),
  testPreset('release-cmake', { inherits: 'default' }),
  testPreset('debug', { inherits: 'default' }),
  testPreset('debug-cmake', { inherits: 'default' }),
  testPreset('dbg', { execution: { jobs: 2 } }),
  testPreset('doxygen', {}),
] +
// Build type + cpp-std + cmake
[testPreset('release-cpp' + std + '-cmake', { inherits: 'default' }) for std in cppStds] +
// Build type + compiler + cmake
[testPreset('release-' + c + '-cmake', { inherits: 'default' }) for c in compilers] +
[testPreset('debug-' + c + '-cmake', { inherits: ['default'] }) for c in compilers] +
// Build type + compiler + cpp-std + cmake (matrix)
std.flattenArrays([
  [testPreset('release-' + c + '-cpp' + std + '-cmake', { inherits: 'default' }) for std in cppStds]
  for c in compilers
]) +
// Spack presets
[testPreset('release-clang-spack', { inherits: 'default' })] +
[testPreset(spackPresetName(env), { inherits: 'default' }) for env in spackEnvironments] +
[testPreset('release-clang-spack-macosx', { inherits: 'default' })] +
[testPreset('debug-clang-spack', { inherits: 'default' })] +
[testPreset('debug-clang-spack-macosx', { inherits: 'default' })] +
[testPreset('release-clang-cpp20-spack', { inherits: 'default' })] +
// Component presets (inherit from default, use 4 jobs, retry failed tests 3 times)
[
  if comp == 'mor' then
    testPresetWithRetry(comp, {
      inherits: 'default',
      execution+: { jobs: 4 },
      filter: {
        include: {
          name: '^(feelpp_test_.*|feelpp_mor_python-tests-(mpi-)?(mor|nirb))$',
        },
      },
    })
  else
    testPresetWithRetry(comp, { inherits: 'default', execution+: { jobs: 4 } })
  for comp in components
] +
// Special presets
[
  testPreset('feelpp-usrlocal', { inherits: 'feelpp' }),
  testPreset('feelpp-cpp20-spack', {}),
  testPreset('feelpp+specx', {}),
  testPreset('mor-dbg', {}),
  testPreset('research-mor', {}),
  testPreset('feelpp-python', {}),
  testPreset('feelpp-python-dbg', {}),
  testPreset('mor_python', {}),
] +
pythonCiTestPresets;

// ============================================================================
// Final Output
// ============================================================================

{
  version: 6,
  cmakeMinimumRequired: {
    major: 3,
    minor: 25,
    patch: 0,
  },
  configurePresets: configurePresets,
  buildPresets: buildPresets,
  testPresets: testPresets,
  workflowPresets: workflowPresets,
  vendor: {
    'example.com/ExampleIDE/1.0': {
      autoFormat: false,
    },
  },
}
