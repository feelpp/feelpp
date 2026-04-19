# Packaging Metadata

`packaging/` is the home of Feel++ packaging metadata and packaging assets.

Current scope:

- Debian and Ubuntu source package trees under `packaging/debian`
- component and distro manifest data under `packaging/manifest`
- `pbuilder` configuration and hooks under `packaging/pbuilder`
- future package-manager metadata such as:
  - `spack/`
  - `guix/`
  - `homebrew/`

This subtree should not contain the Python implementation of the packaging CLI.
That code now belongs under `ops/`.

Current Python tooling:

- implementation: `ops/src/feelpp/pkg`
- preferred CLI: `fpp-pkg`
- compatibility CLI: `feelpp-pkg`

This subtree should stay source-only. Generated artifacts such as:

- `.pytest_cache`
- `__pycache__`
- `*.egg-info`

must not be kept here.
