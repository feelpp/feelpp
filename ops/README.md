# `feelpp-ops`

`ops/` is the home of Feel++ Python operational tooling.

Current commands:

- `fpp-pkg` (preferred)
- `fpp-version` (preferred)

Current Python namespaces:

- `feelpp.pkg`
- `feelpp.ops.common`

Scaffolded sibling namespaces:

- `feelpp.simulate`

Planned later namespaces:

- `feelpp.dataset`

Bootstrap:

```bash
python3 -m venv .venv-ops
. .venv-ops/bin/activate
pip install -e ops
fpp-pkg --help
fpp-version --help
```

Run the `ops` test suite with `pytest`:

```bash
pytest -q ops/tests
```

`fpp-version` uses the repository root [`feelpp.version.cmake`](../feelpp.version.cmake)
as the upstream source of truth and per-component `package_revision` metadata from
[`packaging/manifest/components.toml`](../packaging/manifest/components.toml) to
derive package versions across all configured distros/flavors. Optional
`package_revision_by_dist` overrides let you bump selected distros without
changing the default revision fallback. Use:

```bash
fpp-version show
fpp-version sync
fpp-version sync --dist noble --dry-run
fpp-version revision bump --dist noble --dist trixie --dry-run
fpp-version publications --pretty
fpp-version publications --pretty --rows 10 --since 2026-01-01
fpp-version release 0.111.0-preview.13 --dist noble --dist trixie --dry-run
fpp-version release 0.111.0-preview.13 --dist noble --dist trixie --dry-run --pretty
fpp-version release 0.111.0-preview.13 --dist noble --dist trixie --dry-run --pretty --publications-rows 10 --publications-since 2026-01-01
```

to inspect the desired package versions and then align Debian changelog heads with
that central version state. Upstream GitHub releases stay tied to the upstream
semantic version, while `--dist` lets release validation and notes cover only the
distros shipped in that release. `release --dry-run --pretty` shows the explicit
Feel++ packaging availability notes, install/pull commands for the released
distros, plus GitHub-generated release notes preview, without the lower-level
validation detail blocks. Flavor and distro release versions are sourced from
the packaging catalog in [`.github/plan-ci.json`](../.github/plan-ci.json).
`publications --pretty` harvests recent HAL records from the `FEEL` and
`CEMOSIS` collections and renders them as a Markdown section that can be added
to release notes or discussion posts. Use `--rows` to change how many
publications are shown and `--since` to limit the slice to a release-specific
window. The release command mirrors those controls with
`--publications-rows` and `--publications-since`. The default HAL publication
settings live in [`ops/pyproject.toml`](./pyproject.toml).

The `feelpp.ops.common` namespace is the shared home for cross-tool support
code such as naming, future logging helpers, and execution/runtime helpers.

The Python tooling in `ops/` should stay source-only. Generated artifacts such
as `.pytest_cache`, `__pycache__`, and `*.egg-info` must not be kept here.
