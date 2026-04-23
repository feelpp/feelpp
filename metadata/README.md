# Public IP/APP Metadata

This directory is the public source of truth for Feel++ IP/APP metadata that can be stored in the public `feelpp/feelpp` repository.

The Word DOCX expected by SATT Conectus is not the source of truth here. The canonical public metadata is [`software.public.yml`](./software.public.yml), and public release bundles are generated from that YAML plus repository release metadata.

Private inventor, HR, personal address, payroll, signature, ownership allocation, and confidential company notes must not be stored in this repository. Final APP and SATT Conectus dossier assembly is delegated to the private `cemosis/software-ip` repository.

## Schema

`software.public.yml` uses `schema_version: 1` and requires these top-level sections:

- `software`: public name, short name, description, homepage, and source paths.
- `repository`: public repository identity, default branch, URL, and issue tracker.
- `languages`: public programming language entries with source paths.
- `build_tools`: public build and packaging tools with source paths.
- `public_license`: public license expressions, files, and notes.
- `architecture`: public architecture overview and component path summaries.
- `public_funding`: public funding/grant blocks already present in public metadata.
- `public_outputs`: public citation, CodeMeta, and Zenodo metadata references.
- `metrics`: best-effort source metrics updated by `fpp-ip stats --write-metadata`.
- `private_boundary`: explicit reminder that private dossier assembly is outside this repository.

Unknown or unsafe-to-infer values should be left as `null` rather than guessed.

## Commands

From an editable `ops` install:

```bash
fpp-ip validate-public
fpp-ip show-public
fpp-ip show-public --format json
fpp-ip export-public --format yaml --out metadata/exports/feelpp.public-app.yml
fpp-ip export-public --format json --out metadata/exports/feelpp.public-app.json
fpp-ip stats
fpp-ip stats --write-metadata
```

The committed example export under `metadata/exports/` is generated from public metadata only.
