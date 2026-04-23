# Docker Image Metadata

`packaging/docker/` is the repository-owned home for source-build OCI image
metadata used by `fpp-pkg image ...`.

Current scope:

- Ubuntu and Debian environment metadata under `config/`
- shared environment assets under `assets/`
- base environment templates plus component/full Dockerfile templates under
  `templates/`

Design rules:

- keep only source-of-truth inputs here
- do not commit generated Docker contexts here
- do not mirror the full legacy external Docker repository here
- keep Ubuntu and Debian on the `apt` / component-image path
- keep Spack image generation under `packaging/spack/` plus `ops/`

Current supported variant:

- `feelpp-env`
  - base development image used by the generated Ubuntu/Debian component bake
    graph

Generated bake contexts belong under the packaging job root and are produced by
`fpp-pkg image bake --target ...`.

