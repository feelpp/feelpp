local archs = [
  'native',
  'x86_64',
  'aarch64',
];

local types = [
  'default',
  'asan',
  'spack',
];

local compilers = [
    'clang',
    'gcc',
];
local cpps = [
    'cpp17',
    'cpp20',
    'cpp23',
];
local gpus = [
  'cpu',
  'rocm',
  'cuda',
];

local configs = [
  'Debug',
  'Release',
];

local components = [
    'feelpp',
    'feelpp-core',
    'feelpp-testsuite',
    'feelpp-toolboxes',
    'feelpp-mor',
    'feelpp-python',
];
local cp_generator(component, compiler, cpp, type, gpu, config) =
  {
    name: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
    displayName: component + ' |' + compiler + '|' + cpp + '|' + type + '|' + gpu + '|' + std.asciiLower(config),
    inherits: [cpp, compiler, gpu, type, config,component],
  };

local bp_generator(component, compiler, cpp, type, gpu, config) =
  {
    name: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
    displayName: component + ' |' + compiler + '|' + cpp + '|' + type + '|' + gpu + '|' + std.asciiLower(config),
    configurePreset: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
    configuration: config,
    inherits: "default"
  };

local tp_generator(component, compiler, cpp, type, gpu, config) =
  {
    name: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
    configurePreset: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
    output: { outputOnFailure: true },
    execution: { noTestsAction: 'error', stopOnFailure: true },
  };

local pp_generator(component, compiler, cpp, type, gpu, config) =
  {
    name: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
    steps: [
      {
        type: 'configure',
        name: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
      },
      {
        type: 'build',
        name: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
      },
      {
        type: 'test',
        name: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
      },
      {
        type: 'package',
        name: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
      },
    ],
  };

local wp_generator(component, compiler, cpp, type, gpu, config) =
  {
    name: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
    configurePreset: component + '-' + compiler + '-' + cpp + '-' + type + '-' + gpu + '-' + std.asciiLower(config),
    generators: [
      'TGZ',
    ],
  };

{
  version: 3,
  cmakeMinimumRequired: {
    major: 3,
    minor: 22,
    patch: 0,
  },
  configurePresets: [
    {
      name: 'default',
      hidden: true,
      displayName: 'Default Config',
      description: 'Default build using Ninja Multi-Config generator',
      //generator: 'Ninja Multi-Config',
      generator: "Unix Makefiles",
      binaryDir: "${sourceDir}/build/${presetName}$env{DISTRIBUTION}$env{ARCH}",
      cacheVariables: {
        CMAKE_INSTALL_PREFIX: "${sourceDir}/install/${presetName}$env{DISTRIBUTION}$env{ARCH}",
        CMAKE_EXPORT_COMPILE_COMMANDS: 'ON',
        CMAKE_VERBOSE_MAKEFILE: 'ON',
        FEELPP_ENABLE_VTK: "OFF",
        FEELPP_ENABLE_OPENTURNS: "OFF",
        FEELPP_ENABLE_FMILIB: "OFF",
      },
      environment: {
        NINJA_STATUS: "[run %r|beg %s|fin %f|tot %t|rate %o|time %e]:"
      }
    },
    {
      name: 'asan',
      inherits: 'default',
      hidden: true,
      cacheVariables: {
        CMAKE_CXX_FLAGS_SANITIZE: '-U_FORTIFY_SOURCE -O2 -g -fsanitize=address,undefined -fno-omit-frame-pointer -fno-common',
      },
    },
    {
        name: "spack",
        hidden: true,
        displayName: "spack package manager",
        description: "spack config",
        inherits: [
            "default"
        ],
        cacheVariables: {
            CMAKE_INSTALL_RPATH_USE_LINK_PATH: "ON",
            FEELPP_ENABLE_ANN: "OFF",
            FEELPP_USE_EXTERNAL_CLN: "ON",
            FEELPP_USE_EXTERNAL_GFLAGS: "ON",
            FEELPP_USE_EXTERNAL_GLOG: "ON",
            FEELPP_ENABLE_VTK: "OFF",
            USE_VTK: "OFF",
            FEELPP_ENABLE_OPENTURNS: "OFF"
        },
        environment: {
            VERBOSE: "1"
        }
    },
    {
        name: "cpp17",
        hidden: true,
        description: "Enable compiler C++17",
        cacheVariables: {
            "FEELPP_STD_CPP": "17",
        }
    },
    {
        name: "cpp20",
        hidden: true,
        description: "Enable compiler C++20",
        cacheVariables: {
            "FEELPP_STD_CPP": "20",
        }
    },
    {
        name: "cpp23",
        hidden: true,
        description: "Enable compiler C++23",
        cacheVariables: {
            "FEELPP_STD_CPP": "23",
        }
    },
    {
        name: "gcc",
        hidden: true,
        cacheVariables: {
            CMAKE_C_COMPILER: "gcc",
            CMAKE_CXX_COMPILER: "g++",
        }
    },
    {
        name: "clang",
        hidden: true,
        cacheVariables: {
            CMAKE_C_COMPILER: "clang",
            CMAKE_CXX_COMPILER: "clang++",
        }
    },
    {
        name: 'cpu',
        hidden: true,
        cacheVariables: {
            FEELPP_ENABLE_KOKKOS: 'ON',
            FEELPP_ENABLE_ROCM: 'OFF',
            FEELPP_ENABLE_CUDA: 'OFF',
        /* Additional CPU-specific settings */
        },
    },
    {
        name: 'rocm',
        hidden: true,
        inherits: [
            'feelpp-core-tests-only',
        ],
        cacheVariables: {
            FEELPP_ENABLE_KOKKOS: 'ON',
            FEELPP_ENABLE_ROCM: 'ON',
            FEELPP_ENABLE_CUDA: 'OFF',
            /* Additional ROCm-specific settings */
        },
    },
    {
        name: 'cuda',
        hidden: true,
        inherits: [
            'feelpp-core-tests-only',
        ],
        cacheVariables: {
            FEELPP_ENABLE_ROCM: 'OFF',
            FEELPP_ENABLE_CUDA: 'ON',
            /* Additional CUDA-specific settings */
        },
    },
    {
        name: "Release",
        hidden: true,
        displayName: "Release",
        description: "Plain release build without dependency handling.",
        inherits: [
            "default"
        ],
        cacheVariables: {
            CMAKE_BUILD_TYPE: "Release"
        }
    },
    {
        name: "Debug",
        hidden: true,
        displayName: "Debug|no package manager",
        description: "Plain debug build without dependency handling.",
        inherits: [
            "default"
        ],
        cacheVariables: {
            CMAKE_BUILD_TYPE: "Debug",
            CMAKE_VERBOSE_MAKEFILE: "TRUE",
            CMAKE_MESSAGE_LOG_LEVEL: "VERBOSE",
            H5PP_ENABLE_TESTS: "TRUE",
            H5PP_BUILD_EXAMPLES: "TRUE",
            H5PP_ENABLE_ASAN: "TRUE",
            H5PP_ENABLE_PCH: "FALSE",
            H5PP_ENABLE_CCACHE: "FALSE",
            CMAKE_INTERPROCEDURAL_OPTIMIZATION: "FALSE",
            CMAKE_COMPILE_WARNING_AS_ERROR: "FALSE",
            FEELPP_ENABLE_TOOLBOXES: "ON",
            FEELPP_ENABLE_MOR: "OFF",
            FEELPP_ENABLE_FEELPP_PYTHON: "OFF",
            FEELPP_INSTANTIATION_ORDER_MAX: "1",
            FEELPP_MESH_MAX_ORDER: "1"
        },
        environment: {
            VERBOSE: "1"
        }
    },
    {
        name: "feelpp",
        inherits: [
            "clang",
            "cpp17",
            "Release",
            "default",
        ],
        displayName: "feelpp : clang|cpp17|release",
        description: "Build only the Feel++ library Component",
        cacheVariables: {
            FEELPP_ENABLE_MOR: "OM",
            FEELPP_ENABLE_TOOLBOXES: "ON",
            FEELPP_ENABLE_FEELPP_PYTHON: "ON",
            FEELPP_ENABLE_TESTS: "ON",
            FEELPP_ENABLE_FMILIB: "OFF",
            FEELPP_ENABLE_BENCHMARKS: "OFF",
        }
    },
    {
        name: "feelpp-core-tests-only",
        hidden: true,
        displayName: "feelpp-core only",
        description: "Build only the Feel++ library Component",
        cacheVariables: {
            FEELPP_ENABLE_MOR: "OFF",
            FEELPP_ENABLE_TOOLBOXES: "OFF",
            FEELPP_ENABLE_FEELPP_PYTHON: "OFF",
            FEELPP_ENABLE_TESTS: "ON",
            FEELPP_ENABLE_BENCHMARKS: "OFF",
        }
    },
    {
        name: "feelpp-core",
        inherits: [
            "clang",
            "cpp17",
            "Release",
        ],
        displayName: "feelpp-core : clang|cpp17|release",
        description: "Build only the Feel++ library Component",
        cacheVariables: {
            FEELPP_ENABLE_MOR: "OFF",
            FEELPP_ENABLE_TOOLBOXES: "OFF",
            FEELPP_ENABLE_FEELPP_PYTHON: "OFF",
            FEELPP_ENABLE_TESTS: "OFF",
            FEELPP_ENABLE_BENCHMARKS: "OFF",
        }
    },
    {
        name: "feelpp-testsuite",
        inherits: [
            "clang",
            "cpp17",
            "Release",
        ],
        displayName: "feelpp-testsuite : clang|cpp17|release",
        description: "Build only the Feel++ library Component",
        cacheVariables: {
            FEELPP_COMPONENT: "testsuite",
        }
    },
    {
        name: "feelpp-toolboxes",
        inherits: [
            "clang",
            "cpp17",
            "Release",
        ],
        displayName: "feelpp-toolboxes : clang|cpp17|release",
        description: "Build only the Feel++ Toolboxes Component",
        cacheVariables: {
            FEELPP_COMPONENT: "toolboxes"
        }
    },
    {
        name: "feelpp-mor",
        inherits: [
            "clang",
            "cpp17",
            "Release",
        ],
        displayName: "feelpp-mor : clang|cpp17|release",
        description: "Build only the Feel++ MOR Component",
        cacheVariables: {
            FEELPP_COMPONENT: "mor",
            FEELPP_ENABLE_RESEARCH: "OFF",
            FEELPP_ENABLE_OPENTURNS: "ON"
        }
    },
    {
        name: "feelpp-python",
        inherits: [
            "clang",
            "cpp17",
            "Release",
        ],
        displayName: "feelpp-python : clang|cpp17|release",
        description: "Build only the Feel++ python Component",
        cacheVariables: {
            FEELPP_COMPONENT: "python"
        }
    },
  ] + [cp_generator(component, compiler, cpp, type, gpu, config) for component in components for compiler in compilers for cpp in cpps for type in types for gpu in gpus for config in configs],

buildPresets: [
  {
      name: "default",
      configurePreset: "feelpp",
      jobs: 25
  },
  {
      name: "feelpp",
      configurePreset: "feelpp",
      inherits: "default"
  },
  {
      name: "feelpp-core",
      configurePreset: "feelpp-core",
      inherits: "default"
  },
  {
      name: "feelpp-testsuite",
      configurePreset: "feelpp-testsuite",
      inherits: "default"
  },
  {
      name: "feelpp-toolboxes",
      configurePreset: "feelpp-toolboxes",
      inherits: "default"
  },
  {
      name: "feelpp-mor",
      configurePreset: "feelpp-mor",
      inherits: "default"
  },
  {
      name: "feelpp-python",
      configurePreset: "feelpp-python",
      inherits: "default"
  },
] + [bp_generator(component, compiler, cpp, type, gpu, config) for component in components for compiler in compilers for cpp in cpps for type in types for gpu in gpus for config in configs],
testPresets: [
    {
      name: "default",
      configurePreset: "default",
      filter:{
          include:{
              name:"feelpp"
          }
      },
      output: {
          outputOnFailure: true
      },
      execution: {
          jobs:2,
          noTestsAction: "error",
          stopOnFailure: false,
          repeat: {
              mode: "until-pass",
              count: 3
          }
      }
  },
  {
      name: "feelpp",
      configurePreset: "feelpp",
      inherits:"default",
      execution: {
          timeout: 700
      }
  },
  {
      name: "feelpp-core",
      configurePreset: "feelpp-core",
      inherits:"default",
      execution: {
          timeout: 700
      }
  },
  {
      name: "feelpp-toolboxes",
      configurePreset: "feelpp-toolboxes",
      inherits: "feelpp",
      execution: {
          timeout: 460
      }
  },
  {
      name: "feelpp-mor",
      configurePreset: "feelpp-mor",
      inherits: "feelpp",
      execution: {
          timeout: 240
      }
  }
  {
      name: "feelpp-python",
      configurePreset: "feelpp-python",
      inherits: "feelpp",
      execution: {
          timeout: 700
      }
  },
  {
      name: "feelpp-testsuite",
      configurePreset: "feelpp-testsuite",
      inherits: "feelpp"
  },
] + [tp_generator(component, compiler, cpp, type, gpu, config) for component in components for compiler in compilers for cpp in cpps for type in types for gpu in gpus for config in configs],
//packagePresets: [] + [pp_generator(component, compiler, cpp, type, config) for component in components for compiler in compilers for cpp in cpps for type in types for gpu in gpus for config in configs],
//workflowPresets: [] + [wp_generator(component, compiler, cpp, type, config) for component in components for compiler in compilers for cpp in cpps for type in types for gpu in gpus for config in configs],
}
