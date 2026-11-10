# feeds

Versioned source packages for the [CHD atlas](https://chd-atlas.github.io/).

Each folder under `sources/` reads one upstream source into the atlas's table
format and publishes it as a GitHub release `<source_id>/<version>`: a
`manifest.json` and one gzipped TSV per table. The atlas pins a package in its
`mirrors/sources.yaml` and ingests it with `chd-atlas ingest`; the atlas is the
only thing that can refuse data.

## Setup

```bash
uv sync --python 3.12
```

On a CPU without AVX2 the stock polars wheel dies with "Illegal instruction".
After each `uv sync`, install `polars-runtime-compat==<the polars version in uv.lock>`
into `.venv` and run the gates through `.venv/bin/python -m ...` (a later `uv sync`
or `uv run` removes it again). CI runners have AVX2 and use `uv run`.

## Gates

```bash
ruff check .
ruff format --check .
mypy --strict src
pytest
```

No test and no CI step uses the network.

## The atlas's schemas

`schemas/atlas/` holds the atlas's exported table and manifest JSON Schemas,
read from a git commit (never a working tree); `schemas/atlas/ATLAS_COMMIT`
names it. Refresh with `scripts/sync_atlas_schemas.py <atlas_checkout> [--commit SHA]`.
`feeds.package.write_package` checks every table and the manifest against them
before writing a byte.

## Making a release

1. Merge the recipe to `dev` (the release guard refuses a `recipe_commit` that
   is not on `origin/dev`).
2. Build the package into an empty `dist/`.
3. Dry run: `python -c "from pathlib import Path; from feeds.release import publish; print(publish(Path('dist')))"`.
4. The owner publishes with `dry_run=False`. A tag that exists with identical
   assets is `unchanged`; one with different assets is refused.
