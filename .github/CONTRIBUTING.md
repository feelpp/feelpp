# Contributing to Feel++

Thank you for improving Feel++! This guide aligns with the refreshed coding rules, tooling, and automation so that
all contributions remain consistent across C++20, Python, and HPC workloads.

## License Agreement

By contributing, you agree to license your work under the LGPL or GPL, matching the project license. Only submit
content you are authorized to share.

## Code of Conduct

We follow the [Feel++ Code of Conduct](../CODE_OF_CONDUCT.md). Treat fellow contributors with respect and kindness.

## Getting Started

1. Fork the repository and create a focused branch.
2. Configure the project using the provided CMake preset:
   ```bash
   cmake --preset default -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
   ```
3. Build only what you need:
   ```bash
   cmake --build build/default -j
   ```
4. Enable the Python bindings when relevant: `cmake --build build/default -j --target _core`.
5. Run focused tests (`ctest -R <regex>` or `pytest` in `python/pyfeelpp`).

## Formatting & Style

Feel++ enforces a single formatting profile stored in `.clang-format` (Allman braces, 4-space indentation,
column limit 100, pointers aligned with the type). Before committing:

```bash
clang-format -i path/to/file.cpp
```

### Exclusions

Do **not** format or tidy the vendor trees: `third_party/` and `external/`. The pre-commit hooks and CI configuration
already skip those paths. If you need to exclude additional generated code, add patterns to `.pre-commit-config.yaml`.

## Static Analysis

Clang-Tidy is configured via `.clang-tidy` to run targeted checks:
- `bugprone-*`, `performance-*`, `modernize-*` (treated as errors)
- `readability-identifier-naming` enforcing the Feel++ naming conventions

Run clang-tidy before sending a PR:

```bash
cmake --build build/default -j            # ensure the project is up-to-date
ninja -C build/default clang-tidy        # if a tidy target exists
# or fallback
run-clang-tidy -p build/default path/to/file.cpp
```

If your generator does not emit a tidy target, create one in your local `CMakeLists.txt` using `clang-tidy` tooling.

## Pre-Commit Hooks

Install and run the bundled hooks for consistent results:

```bash
pip install pre-commit
pre-commit install
pre-commit run --all-files
```

Hooks executed:
- `clang-format` (skips `third_party/` and `external/`)
- `codespell` for lightweight spelling checks

Expect pre-commit to block commits when formatting or spelling issues are present. Fix the reported files and rerun.

## Continuous Integration

GitHub Actions validates formatting and static analysis on every PR:
- **clang-format-check**: ensures no diff is introduced by clang-format.
- **clang-tidy**: runs targeted diagnostics using `compile_commands.json`.

Your PR must pass these jobs in addition to any existing build/test pipelines.

## Commit & PR Etiquette

- Use imperative commit subjects: `component: concise summary`.
- Keep commits focused; avoid mixing tooling churn with feature changes.
- Reference issues with `Fixes #123` or `Refs #123` in the PR description when applicable.
- In the PR body, describe the testing performed and mention any skipped checks or TODOs.

## Reporting Issues

File issues on GitHub with clear reproduction steps, compiler/tool versions, and minimal Feel++ examples.
Include build logs and relevant snippets; a failing test case earns priority.

## Need Help?

Reach out via GitHub Discussions or the Feel++ community channels. Provide context (OS, compiler, steps taken) so we can
assist efficiently.

