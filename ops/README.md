# `feelpp-ops`

`ops/` is the home of Feel++ Python operational tooling.

Current commands:

- `fpp-pkg` (preferred)

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
```

Run the package-tool test suite with `pytest`:

```bash
pytest -q ops/tests/pkg
```

The `feelpp.ops.common` namespace is the shared home for cross-tool support
code such as naming, future logging helpers, and execution/runtime helpers.

The Python tooling in `ops/` should stay source-only. Generated artifacts such
as `.pytest_cache`, `__pycache__`, and `*.egg-info` must not be kept here.
