# Feel++ I/O and Logging: MPI-Safe Design

## Overview

Feel++ provides MPI-aware I/O streams and logging that work safely throughout the entire program lifecycle, including:
- Before MPI initialization
- During MPI execution
- After MPI finalization
- In non-MPI programs

## I/O Streams: `Feel::cout`, `Feel::cerr`, `Feel::clog`

### Design

The `MasterStream` class (defined in `feelcore/feelio.hpp`) wraps standard C++ streams with MPI-aware behavior:

```cpp
Feel::cout << "Message" << std::endl;  // Only outputs on master rank when MPI is active
Feel::cerr << "Error" << std::endl;    // Only outputs on master rank when MPI is active
```

### Behavior Based on MPI State

| MPI State | Behavior | Rationale |
|-----------|----------|-----------|
| **Not initialized** | Outputs on all processes | MPI not available yet, safe to output |
| **Initialized & Active** | Outputs only on master rank | Avoids duplicate output from all ranks |
| **Finalized** | Outputs on all processes | Cannot safely call MPI functions, avoid crash |
| **No MPI** | Always outputs | Non-MPI programs work normally |

### Implementation Details

The `shouldOutput()` method checks MPI state without causing crashes:

```cpp
bool shouldOutput() const noexcept
{
    // Safe checks that never call MPI after finalization
    if ( !Environment::initialized() || Environment::finalized() )
        return true;  // Always output when MPI not active
    
    // MPI active: check if master rank
    return wc->isMasterRank();
}
```

### Usage Examples

```cpp
// Example 1: Before Environment
int main(int argc, char** argv)
{
    // Works fine - outputs on all processes
    Feel::cout << "Starting application" << std::endl;
    
    Environment env(argc, argv);
    
    // Only master rank outputs
    Feel::cout << "Environment initialized" << std::endl;
    
    return 0;  // After env destroyed, Feel::cout still works
}

// Example 2: After Environment destruction
int main(int argc, char** argv)
{
    {
        Environment env(argc, argv);
        Feel::cout << "Inside scope" << std::endl;
    }
    // MPI finalized, but Feel::cout still safe
    Feel::cout << "After Environment" << std::endl;
    return 0;
}
```

## Logging: Google glog Integration

### System Overview

Feel++ uses [Google's glog library](https://github.com/google/glog) for logging:

```cpp
LOG(INFO) << "Informational message";
LOG(WARNING) << "Warning message";
LOG(ERROR) << "Error message";
VLOG(1) << "Verbose message at level 1";
VLOG(2) << "Verbose message at level 2";
```

### Key Functions

#### `Environment::startLogging()`

Initializes the logging system:

```cpp
Environment::startLogging();
LOG(INFO) << "This goes to log file";
```

#### `Environment::stopLogging(bool remove = false)`

Safely shuts down logging with MPI-aware cleanup:

```cpp
// Just stop logging
Environment::stopLogging();

// Stop logging and delete log files (master rank only if MPI active)
Environment::stopLogging(true);
```

**MPI Safety**: After MPI finalization (e.g., after `PetscFinalize()`), `stopLogging()` will run cleanup on all ranks to avoid calling `isMasterRank()` which would crash.

### Logging Lifecycle

```
Program Start
    ↓
[Optional: Early logging to stderr]
    ↓
Environment Construction
    ↓
startLogging() → Logs to files in logsRepository()
    ↓
Application Execution (LOG, VLOG work normally)
    ↓
PetscFinalize() → MPI_Finalize called
    ↓
stopLogging() → Safe cleanup (no MPI calls)
    ↓
Environment Destruction
    ↓
Program End (Feel::cout still works)
```

### Verbosity Control

Set verbosity level for VLOG macros:

```cpp
// Command line
./app --v=2  // Show VLOG(1) and VLOG(2)

// Programmatically
Environment::setVerboseLevel(2);
```

### Log File Location

Logs are stored in the application repository:

```cpp
std::string log_path = Environment::logsRepository();
// Typically: $FEELPP_DIR/app_name/np_N/logs/
```

## Common Patterns

### Pattern 1: Standard Application

```cpp
int main(int argc, char** argv)
{
    Environment env(argc, argv, makeAbout(), makeOptions());
    
    // I/O and logging work normally
    Feel::cout << "Application started" << std::endl;
    LOG(INFO) << "Logging initialized";
    
    // Your application code
    
    return 0;  // Automatic cleanup, no MPI errors
}
```

### Pattern 2: Early Messages

```cpp
int main(int argc, char** argv)
{
    // Before Environment: outputs on all processes
    Feel::cout << "Pre-initialization" << std::endl;
    
    Environment env(argc, argv);
    
    // After Environment: master rank only
    Feel::cout << "Initialization complete" << std::endl;
    
    return 0;
}
```

### Pattern 3: Manual Log Control

```cpp
int main(int argc, char** argv)
{
    Environment env(argc, argv);
    
    Environment::startLogging();
    LOG(INFO) << "Starting computation";
    
    // ... computation ...
    
    Environment::stopLogging();
    LOG(INFO) << "This won't be logged";
    
    return 0;
}
```

## Debugging Tips

### Enable Verbose Logging

```bash
# Show all VLOG up to level 3
./app --v=3

# Show specific modules
./app --vmodule=environment=2,repository=1
```

### Log to stderr

```bash
# Don't create log files, output to terminal
export GLOG_logtostderr=1
./app
```

### MPI Debugging

If you see "MPI_Comm_rank() called after MPI_FINALIZE":

1. Check if `Feel::cout`/`Feel::cerr` is used after Environment destruction
2. Check if `isMasterRank()` is called directly after MPI finalization
3. Verify `stopLogging()` is called before MPI finalization if possible

## Technical Implementation

### Thread Safety

`MasterStream::shouldOutput()` is marked `noexcept` and uses exception handling to safely recover from MPI errors during destruction.

### Exception Safety

All MPI calls in I/O and logging infrastructure use try-catch blocks to prevent crashes during program cleanup.

### Backward Compatibility

The redesigned `MasterStream` maintains 100% backward compatibility:
- Existing code continues to work
- Constructor signature unchanged
- All public methods preserved
- Behavior identical during normal execution

The only difference is improved safety during program initialization and cleanup phases.

## Migration Guide

### From std::cout to Feel::cout

```cpp
// Before
if (Environment::isMasterRank())
    std::cout << "Message" << std::endl;

// After
Feel::cout << "Message" << std::endl;
```

### From Direct MPI Checks

```cpp
// Before (unsafe after MPI finalization)
if (worldComm()->isMasterRank())
    doSomething();

// After (safe throughout lifecycle)
if (Environment::initialized() && !Environment::finalized() && 
    worldComm()->isMasterRank())
    doSomething();
```

## Performance Notes

- `shouldOutput()` is inline and very fast (few CPU cycles)
- MPI state checks are cached where possible
- No performance impact during normal execution
- Minimal overhead even when output is suppressed

## Testing

See `doc/manual/examples/test_feelio.cpp` for comprehensive tests of all scenarios.

## Summary

Feel++ I/O and logging systems are designed to be:
- **MPI-aware**: Automatically handle master-only vs. all-rank output
- **Lifecycle-safe**: Work before init, during execution, and after finalization
- **Exception-safe**: Recover gracefully from errors
- **Backward-compatible**: No changes needed to existing code
- **Well-documented**: Clear semantics and behavior

This ensures robust, crash-free programs even during complex initialization/cleanup phases with MPI.
