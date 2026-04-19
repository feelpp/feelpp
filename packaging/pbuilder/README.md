This directory contains the static `pbuilder-dist` configuration used by the
Feel++ Debian/Ubuntu packaging jobs.

- `pbuilderrc` defines shared pbuilder behavior.
- `hooks/` contains the chroot customization hooks copied into `/hooks`.
- `hooks/keyrings/*.gpg.b64` contains the external apt repository keys in
  base64-encoded, dearmored form consumed by the hooks.

`feelpp-pkg pbuilder prepare` derives the heavy preload package set directly
from the active in-tree `packaging/debian/*/<dist>/debian/control` files.
`debian/control` is the source of truth. The runtime preload hook is generated
from those source-package build dependencies after filtering out internal
Feel++ binary packages that are provided later by the local chained-build repo.

The seeded-base prepare path is intentionally retryable: if the preload fails in
the middle of a long dependency download, the partial apt archive cache is kept
under the persistent pbuilder state root and the next prepare attempt resumes
from that cache instead of starting from zero.

Dynamic repository selection still comes from the wrapper in
`feelpp/tools/scripts/pkg/feelpp_pkg_common.sh`, which injects the current
channel and distro policy through `OTHERMIRROR`.
