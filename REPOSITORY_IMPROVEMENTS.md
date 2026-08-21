# Feel++ Repository System Improvements

## Overview

This document describes the improvements made to the Feel++ Repository system, focusing on enhancing flexibility and testability of the database/results directory configuration.

## Changes Summary

### 1. New `Location::custom` Type

Added a new location type to support dynamic repository path determination through callback functions.

**File: `feelpp/feel/feelcore/repository.hpp`**

```cpp
enum class Location {
    unknown=0,    ///! unknown directory
    global=10,    ///! global repository
    relative,     ///! relative to current directory
    absolute,     ///! absolute directory given 
    git,          ///! relative to git repository
    custom        ///! custom location determined by a callback function (NEW)
};
```

### 2. Custom Location Callback Support

**File: `feelpp/feel/feelcore/repository.hpp`**

Added callback function support to `Repository::Config`:

```cpp
struct Config {
    // ... existing fields ...
    
    /// Custom callback to determine repository location
    std::function<fs::path()> custom_location_callback;
    
    // New constructor
    Config( fs::path d, std::function<fs::path()> callback );
};
```

**File: `feelpp/feel/feelcore/repository.cpp`**

Implemented callback invocation in `Repository::configure()`:

```cpp
if ( config_.location == Location::custom )
{
    if ( config_.custom_location_callback )
    {
        try {
            root_ = config_.custom_location_callback();
        }
        catch ( const std::exception& e ) {
            // Fallback to default on error
        }
    }
}
```

### 3. Helper Function

**File: `feelpp/feel/feelcore/repository.hpp`**

New helper function for creating custom repositories:

```cpp
inline Repository::Config customRepository( 
    std::string reldir, 
    std::function<fs::path()> callback )
{
    return Repository::Config(fs::path(reldir), callback);
}
```

### 4. Repository Methods

Added `isCustom()` method to check if repository uses custom location:

```cpp
bool isCustom() const { return config_.location == Location::custom; }
```

## Usage Examples

### C++ Usage with Environment

#### Using Custom Repository in Feel++ Applications

The most common use case is to provide a custom repository configuration when initializing the Environment:

```cpp
#include <feel/feelcore/environment.hpp>

int main(int argc, char** argv)
{
    using namespace Feel;
    
    // Define options for your storage configuration
    po::options_description opts("Storage options");
    opts.add_options()
        ("storage.base", po::value<std::string>()->default_value("/tmp"), 
         "base storage directory");
    
    // Create custom repository that uses parsed options
    auto config = customRepository("fallback", []() {
        std::string base = Feel::soption(_name="storage.base");
        std::string user = Feel::findUser();
        return fs::path(base) / user / "feelpp-results";
    });
    
    // Initialize Environment with custom repository
    Environment env(_argc=argc, _argv=argv,
                   _desc=opts,
                   _config=config);
    
    // Repository is now configured and accessible
    std::cout << "Results in: " << Environment::appRepository() << std::endl;
    
    return 0;
}
```

#### Simple Custom Location

```cpp
auto config = customRepository("fallback", []() {
    return fs::path("/tmp/my-custom-repo");
});

Repository repo(config);
repo.configure();
```

#### Dynamic Location Based on Options

```cpp
auto config = customRepository("fallback", []() {
    // Access Feel++ options within the lambda
    std::string basedir = Feel::soption(_name="custom.basedir");
    std::string subdir = Feel::soption(_name="custom.subdir");
    return fs::path(basedir) / subdir;
});

Repository repo(config);
repo.configure();
```

#### Environment Variable Based Location

```cpp
auto config = customRepository("fallback", []() {
    const char* env_path = getenv("FEELPP_CUSTOM_ROOT");
    if (env_path)
        return fs::path(env_path) / "results";
    else
        return fs::path("/tmp/feelpp-fallback");
});
```

### Python Usage

```python
import feelpp.core as fppc

# Define a Python callback
def compute_custom_path():
    base = "/tmp/my-python-repo"
    return pathlib.Path(base) / "subdir"

# Create custom repository config
config = fppc.customRepository("fallback-dir", compute_custom_path)

# Use it
repo = fppc.Repository(config)
repo.configure()

assert repo.isCustom()
print(f"Repository root: {repo.root()}")
```

## Testing

### New Test Suite

**File: `testsuite/feelcore/test_repository.cpp`**

Comprehensive test suite covering:

1. **test_repository_default** - Basic construction
2. **test_repository_global** - Global location
3. **test_repository_relative** - Relative location
4. **test_repository_absolute** - Absolute location
5. **test_repository_git** - Git-based location
6. **test_repository_git_error** - Git error handling
7. **test_repository_custom** - Custom location with options
8. **test_repository_custom_env** - Custom location with environment variables
9. **test_repository_custom_error** - Custom callback error handling
10. **test_repository_structure** - Directory structure and appenders
11. **test_repository_reconfig** - Reconfiguration
12. **test_repository_user_info** - User information
13. **test_repository_cd** - Change directory functionality

### Running Tests

```bash
# Build and run all repository tests
cd build/default
ctest -R test_repository -V

# Or run specific tests
./testsuite/feelcore/feelpp_test_repository --run_test=repository/test_repository_custom
```

### Python Tests

**File: `python/pyfeelpp/tests/test_core.py`**

Added `test_repository_custom()` to demonstrate Python callback usage.

```bash
# Run Python tests
cd python/pyfeelpp
pytest tests/test_core.py::test_repository_custom -v -s
```

## Benefits

### 1. Flexibility

- **Dynamic Configuration**: Repository location can be determined at runtime based on:
  - Command-line options
  - Environment variables
  - Configuration files
  - System state
  - User preferences

### 2. Use Cases

- **Multi-tenant Systems**: Different users/projects can have isolated repositories
- **Cluster Environments**: Repository can be placed on appropriate storage (scratch, home, project)
- **CI/CD Pipelines**: Dynamic paths based on job/branch/commit
- **Testing**: Isolated test environments without conflicts

### 3. Backward Compatibility

All existing location types (`global`, `relative`, `absolute`, `git`) remain unchanged. The new `custom` type is purely additive.

## Error Handling

The custom callback implementation includes robust error handling:

1. **Exception Catching**: If callback throws, repository falls back to safe default
2. **Logging**: Errors are logged for debugging
3. **Graceful Degradation**: System continues to function even if callback fails

## Environment Integration

### Initialization Flow

The custom repository integrates seamlessly with Feel++ Environment initialization:

1. **Environment Constructor**: Repository is created with custom config but NOT configured yet
2. **Options Parsing**: Command-line and config file options are parsed
3. **changeRepository()**: Called automatically after options parsing
4. **Callback Invocation**: Custom location callback is executed with access to all parsed options
5. **Directory Creation**: Repository directories are created based on computed path

```
Environment Construction Flow:
┌─────────────────────────────────────────────────────┐
│ 1. Create Repository with Config                    │
│    - If custom: defer configure()                   │
│    - If non-custom: configure() immediately         │
└────────────────────┬────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────┐
│ 2. Parse Options (doOptions)                        │
│    - Command-line arguments                         │
│    - Config files                                   │
│    - All Feel++ options available                   │
└────────────────────┬────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────┐
│ 3. changeRepository() called                        │
│    ┌──────────────────────────────────────────────┐ │
│    │ If custom location:                          │ │
│    │   • Invoke callback()                        │ │
│    │   • Callback accesses options via soption()  │ │
│    │   • Directory param IGNORED                  │ │
│    │                                              │ │
│    │ If non-custom:                               │ │
│    │   • Use directory parameter                  │ │
│    │   • Standard location logic                  │ │
│    └──────────────────────────────────────────────┘ │
└────────────────────┬────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────┐
│ 4. Repository configured and ready                  │
│    - Directories created                            │
│    - Environment::appRepository() accessible        │
└─────────────────────────────────────────────────────┘
```

This ensures that custom callbacks have full access to parsed options when computing the repository path.

### Key Points

- Custom location callbacks are **deferred** until after options are parsed
- Callbacks can safely use `Feel::soption()`, `Feel::ioption()`, etc.
- Non-custom locations (global, relative, absolute, git) are configured immediately
- The repository can be reconfigured at runtime using `Environment::changeRepository()`
- **Important**: For custom locations, the callback **always takes precedence** over any directory parameter passed to `configure()` or `changeRepository()`

### Custom Location Precedence

When using `Location::custom`, the callback function is the **sole source of truth** for the repository path:

```cpp
// The callback determines the path
auto config = customRepository("fallback", []() {
    return fs::path("/my/computed/path");
});

Environment env(_config=config);

// Even if changeRepository() is called with a directory,
// the callback's result is used for custom locations
Environment::changeRepository(_directory="ignored-for-custom");

// Result: /my/computed/path (from callback)
// NOT: ignored-for-custom (directory parameter ignored)
```

This design ensures that:
1. Custom logic isn't accidentally overridden
2. Options-based path computation remains consistent
3. Directory parameters only affect non-custom locations

## Design Considerations

### Why Lambda/Callback Approach?

1. **Late Binding**: Repository location is determined *after* options are parsed
2. **Flexibility**: User code can implement any logic
3. **Encapsulation**: Repository class doesn't need to know about all possible location strategies
4. **Testability**: Easy to inject different behaviors for testing
5. **Access to Context**: Callbacks execute when full application context (options, environment) is available

### Thread Safety

The callback is invoked during `Repository::configure()`, which should be called once during initialization. Not designed for concurrent reconfiguration.

### Python Integration

Python callbacks are properly wrapped and converted using pybind11, allowing natural Python code to compute paths.

## Future Enhancements

Potential future improvements:

1. **Async Callbacks**: Support for asynchronous path determination
2. **Validation**: Optional path validation before acceptance
3. **Caching**: Cache computed paths for repeated configure calls
4. **Templates**: Pre-defined templates for common patterns
5. **Policy Classes**: Alternative to lambdas for complex location strategies

## Migration Guide

### From Manual Path Construction

**Before:**
```cpp
std::string base = soption("basedir");
std::string sub = soption("subdir");
fs::path full_path = fs::path(base) / sub;
Repository::Config config(full_path, Location::absolute);
```

**After:**
```cpp
auto config = customRepository("fallback", []() {
    return fs::path(soption("basedir")) / soption("subdir");
});
```

### From Environment Variables

**Before:**
```cpp
const char* path = getenv("FEELPP_REPO");
Repository::Config config(
    path ? fs::path(path) : fs::path("/tmp/default"),
    Location::absolute
);
```

**After:**
```cpp
auto config = customRepository("fallback", []() {
    const char* path = getenv("FEELPP_REPO");
    return path ? fs::path(path) : fs::path("/tmp/default");
});
```

## References

- **Header**: `feelpp/feel/feelcore/repository.hpp`
- **Implementation**: `feelpp/feel/feelcore/repository.cpp`
- **Tests**: `testsuite/feelcore/test_repository.cpp`
- **Python Bindings**: `python/pyfeelpp/feelpp/core/core.cpp`
- **Python Tests**: `python/pyfeelpp/tests/test_core.py`

## Authors

- Christophe Prud'homme <christophe.prudhomme@feelpp.org>

## License

LGPL 2.1+ (same as Feel++ library)
