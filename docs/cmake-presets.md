# CMake Presets Generation

This directory uses **Jsonnet** to generate `CMakePresets.json` from a maintainable template.

## Why Jsonnet?

The previous `CMakePresets.json` had **91 presets in 883 lines** with massive repetition. Using Jsonnet:
- ✅ Reduces to ~600 lines of maintainable code
- ✅ Eliminates repetition through functions and loops
- ✅ Easy to add new compiler/C++ standard/component combinations
- ✅ Single source of truth for configuration patterns
- ✅ Automatically generates matrix combinations

## Setup

### Install Jsonnet

**Ubuntu/Debian:**
```bash
sudo apt install jsonnet
```

**macOS:**
```bash
brew install jsonnet
```

**From source (Go):**
```bash
go install github.com/google/go-jsonnet/cmd/jsonnet@latest
```

## Usage

### Generate CMakePresets.json

```bash
./scripts/generate-presets.sh
```

Or manually:
```bash
jsonnet CMakePresets.jsonnet | python3 -m json.tool --indent 4 > CMakePresets.json
```

### Workflow

1. **Edit** `CMakePresets.jsonnet` (the source template)
2. **Generate** with `./scripts/generate-presets.sh`
3. **Test** with `cmake --preset <preset-name>`
4. **Commit** both `.jsonnet` and generated `.json` files

## Structure of CMakePresets.jsonnet

### Configuration Arrays
```jsonnet
local compilers = ['gcc', 'clang'];
local cppStds = ['20', '23'];
local components = ['feelpp', 'testsuite', 'toolboxes', 'mor', 'python'];
```

### Base Presets
- `default` - Base configuration
- Hidden mixins: `cpp20`, `cpp23`, `gcc`, `clang`, `spack`, `usrlocal`, `warnings`

### Generated Presets
1. **Build types**: `release`, `debug`, `dev`
2. **Matrix combinations**: 
   - `{buildType}-{compiler}-cpp{std}-cmake`
   - Example: `release-clang-cpp23-cmake`
3. **Component builds**: `feelpp`, `toolboxes`, `mor`, `python`, `testsuite`
4. **Special variants**: `feelpp-usrlocal`, `mor-dbg`, `feelpp+specx`

### Adding New Combinations

**Add a new compiler:**
```jsonnet
local compilers = ['gcc', 'clang', 'intel'];  // Add 'intel'
```

**Add a new C++ standard:**
```jsonnet
local cppStds = ['20', '23', '26'];  // Add '26'
```

**Add a new component:**
```jsonnet
local components = ['feelpp', 'testsuite', 'toolboxes', 'mor', 'python', 'newcomp'];

local componentCacheVars = {
  // ... existing ...
  newcomp: {
    FEELPP_COMPONENT: 'newcomp',
    FEELPP_ENABLE_NEWCOMP: 'ON',
  },
};
```

The generation script will automatically create all combinations!

## Benefits

| Aspect | Before | After |
|--------|--------|-------|
| Lines of code | 883 | ~600 |
| Number of presets | 91 | Same (auto-generated) |
| Maintenance | Manual copy-paste | Change once, apply everywhere |
| Adding new compiler | ~30 manual edits | Add 1 line to array |
| Adding C++ std | ~40 manual edits | Add 1 line to array |
| Type safety | JSON (no validation) | Jsonnet (validated) |
| Documentation | Comments scattered | Structure IS documentation |

## Git Workflow

The generated `CMakePresets.json` is committed to git so users without jsonnet can still use presets:

```bash
# After editing CMakePresets.jsonnet
./scripts/generate-presets.sh
git add CMakePresets.jsonnet CMakePresets.json
git commit -m "Update CMake presets: add C++26 support"
```

## Customization

See `CMakePresets.jsonnet` for the full template. Key sections:

- **Lines 1-10**: Configuration arrays (compilers, standards, components)
- **Lines 40-90**: Base preset definitions
- **Lines 100-200**: Mixin presets (hidden)
- **Lines 210-350**: Preset generation functions
- **Lines 360-450**: Component-specific configurations
- **Lines 460+**: Matrix generation and output

## Troubleshooting

**Error: jsonnet command not found**
```bash
sudo apt install jsonnet  # or see Install section above
```

**Validation failed**
```bash
jsonnet --lint CMakePresets.jsonnet
```

**Compare before/after**
```bash
# Backup current
cp CMakePresets.json CMakePresets.json.backup
# Generate new
./scripts/generate-presets.sh
# Compare
diff CMakePresets.json.backup CMakePresets.json
```

## Related Documentation

- [CMake Presets Documentation](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html)
- [Jsonnet Language](https://jsonnet.org/)
- [Jsonnet Tutorial](https://jsonnet.org/learning/tutorial.html)
