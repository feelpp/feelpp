# Site Include Templates

Files in this directory are reserved for documenting site-specific override
patterns for Spack configuration.

Examples of data that belong here rather than in shared repository manifests:

- local compiler externals
- site mirrors
- site package preferences
- accelerator architecture policy tied to one machine or cluster

Phase 0 keeps this directory as documentation-only scaffolding.

Phase 1 adds example files that show the kind of configuration that should live
outside the shared repository-owned environments:

- `compilers.example.yaml`
- `packages.example.yaml`
- `mirrors.example.yaml`
