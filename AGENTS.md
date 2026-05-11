Purpose
- Provide a concise, repo-specific playbook for human and ai agents: how to build, test, patch, and validate across C++ core and Python bindings, with specific guidance for RemoteData/Girder/GitHub/CKAN.

Repository Overview
- C++ core: under `feelpp/feel/...` builds a shared lib `feelpp` and related components.
- Python bindings: under `python/pyfeelpp/feelpp/*` producing `_core` and other py modules installed into `feelpp/core`.
- Tests:
  - Python tests: `python/pyfeelpp/tests/*.py` (includes RemoteData tests).
  - C++ tests: under `testsuite/**`.
  - Tools/CLI sources: `feelpp/tools/**` (includes RemoteData CLI).

Build
- Presets: defined in `CMakePresets.json`. Default build directory: `build/default`.
- Configure (if not already configured):
  - `cmake --preset default`
- Build the core lib and Python modules:
  - `cmake --build build/default -j`
  - Python core target (when needed): `cmake --build build/default -j --target _core`
- Environment hints:
  - Some CI/local sandboxes do not resolve passwd entries for the UID; set `HOME` to a writable temp if needed: `export HOME=$(mktemp -d)`.

Python Tests
- From `python/pyfeelpp`:
  - `export PYTHONPATH=$(pwd)/../../build/default/python/pyfeelpp:$PYTHONPATH`
  - `export LD_LIBRARY_PATH=$(pwd)/../../build/default/feelpp/feel:$LD_LIBRARY_PATH`
  - Run all: `pytest -q -v -s`
  - Focused examples:
    - RemoteData Girder safe test: `pytest -q tests/test_remotedata.py -k girder_operations_safe -v -s`
    - RemoteData comprehensive: `pytest -q tests/test_remotedata.py -v -s`
- Diagnostics:
  - Log gflags/glog to stderr: `env GLOG_logtostderr=1 <pytest cmd>`
  - If `getpwuid()` errors: `export HOME=$(mktemp -d)` and rerun.

CLI (Tools) – RemoteData
- Sources under `feelpp/tools/remotedata/`.
- Typical scenarios covered by Python tests have CLI analogs (download, contents, path lookup, uploads).
- For CI safety, prefer unauthenticated operations by default; authenticated operations require env vars below.

Environment & Secrets
- RemoteData integrations (never log/print these):
  - Girder: `FEELPP_GIRDER_API_KEY`, `FEELPP_GIRDER_TOKEN`
  - CKAN: `FEELPP_CKAN_API_KEY`, `FEELPP_CKAN_URL`, `FEELPP_CKAN_ORGANIZATION`
  - GitHub: `FEELPP_GITHUB_TOKEN`
- Progress verbosity control (Python & CLI):
  - Support QUIET/NORMAL/VERBOSE/DEBUG. Default to QUIET in automated testing to avoid token exposure.

Agent Operating Rules
- Scope changes narrowly to the task; match existing style.
- Fix root causes rather than masking failures.
- Use `apply_patch` for edits; do not `git commit` unless requested.
- Read files in ≤250 line chunks; prefer `rg` for code search.
- Use `update_plan` for multi-step work; exactly one step in progress.
- Send a short preamble before grouped tool calls.

Validation Philosophy
- Rebuild only necessary targets (e.g., `feelpp`, `_core`).
- Start with the most focused tests that exercise the changes, then broaden.
- If network is restricted, limit to public/unauthed operations or skip and document.

Troubleshooting
- Abort in RemoteData Girder parsing:
  - Ensure regex parsing doesn’t abort on valid descriptors; prefer safe parsing and graceful errors.
- Missing passwd/home:
  - If you see `getpwuid(): uid not found`, set `HOME` to a writable path (e.g., `/tmp` via `mktemp -d`).
- Token visibility in logs:
  - Use QUIET or ensure debug printing is gated behind DEBUG level; never print full tokens.

Preferred Test Mix (Remotedata)
- Python: keep a few high-signal tests that mirror CLI functionality:
  - GitHub basic download/contents
  - URL download + invalid URL
  - Girder safe contents/path lookup
  - Progress control and token suppression (QUIET/NORMAL/VERBOSE/DEBUG)
- CLI: add 2–4 smoke tests executed in CI that do not require secrets (e.g., public GitHub/URL operations); keep authenticated CLI tests optional behind secrets.

Commit & Branch Policy
- Do not create branches/commits unless explicitly requested by a maintainer.
- Provide a brief change summary with file paths and next steps after significant edits.

