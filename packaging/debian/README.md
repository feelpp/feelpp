# Feel++ Debian Packaging

This directory is the in-tree home for Feel++ Debian/Ubuntu packaging metadata.

Current scope:

- Ubuntu Noble (`noble`)
- Ubuntu Resolute (`resolute`)
- `feelpp`
- `feelpp-toolboxes`
- `feelpp-mor`

The standalone `feelpp-python` source package is intentionally not kept here.
Python bindings belong to their owning source packages instead of being shipped
as a separate packaging tree.

Current Python binary package ownership:

- `feelpp` ships `python3-feelpp`
- `feelpp-toolboxes` ships `python3-feelpp-toolboxes` and the toolbox-specific `python3-feelpp-toolboxes-*` packages that are actually built
- `feelpp-mor` ships `python3-feelpp-mor`

Layout:

- `packaging/debian/feelpp/noble`
- `packaging/debian/feelpp/resolute`
- `packaging/debian/feelpp-toolboxes/noble`
- `packaging/debian/feelpp-toolboxes/resolute`
- `packaging/debian/feelpp-mor/noble`
- `packaging/debian/feelpp-mor/resolute`

These trees were imported from the legacy `feelpp.pkg` repository so packaging changes can now be reviewed in the main Feel++ source repository.

The legacy `feelpp/tools/scripts/pkg/feelpp_pkg.sh` script now reads metadata from this directory by default instead of cloning `feelpp.pkg`.

Additional Debian and Ubuntu distro trees can be imported here incrementally as the packaging refactor proceeds.
