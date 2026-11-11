# feeds

Versioned source packages for the [CHD atlas](https://chd-atlas.github.io/).

Each folder under `sources/` reads one upstream source into the atlas's table
format and publishes it as a GitHub release `<source_id>/<version>`: a
`manifest.json` and one gzipped TSV per table. The atlas pins a package in its
`mirrors/sources.yaml` and ingests it with `chd-atlas ingest`; the atlas's verdict
is the one that counts.

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
mypy --strict src sources
pytest
```

No test uses the network (the suite passes in a network-less namespace); CI only
installs the locked dependencies.

## The atlas's schemas

`schemas/atlas/` holds the atlas's exported table and manifest JSON Schemas,
read from a git commit (never a working tree); `schemas/atlas/ATLAS_COMMIT`
names it. Refresh with `scripts/sync_atlas_schemas.py <atlas_checkout> [--commit SHA]`.
`feeds.package.write_package` checks every table and the manifest against them
before writing a byte.

## Conventions

Every recipe fills the same manifest fields the same way.

- **`inputs[].url`** is where the hashed bytes can be fetched again, when they can.
  A cluster intermediate that is not published (Farah's pseudobulk, HGNC's hvantk
  export) gets a permalink to the code that made it instead.
- **`inputs[].retrieved_on`** is the day the bytes were obtained or made: the
  operator's download day, the date in hvantk's sidecar, or the cluster job's date.
- **`built_by.recipe_commit`** is read from git by the build, never typed. The
  build refuses a checkout with uncommitted changes under `src`, `sources`,
  `schemas`, `pyproject.toml` or `uv.lock`, and refuses a caller-supplied value
  that differs. **`built_by.host`** is the machine that ran the step the operator
  names: the cluster node for a job, the workstation for a recipe-only build.
  **`built_by.hvantk_release`** (hvantk-based recipes) is the full hvantk commit
  that ran.
- **The revision rule.** A published version is permanent: one version, one set of
  bytes (F13). `revision` is a required argument. It is 1 for the first build of an
  input and 2, 3, ... (`-r2`, ...) for any later build whose `manifest.json` bytes
  differ. That includes a recipe change, and also any rebuild that changes
  `built_by` (a later `recipe_commit`, another host or job) or `retrieved_on`,
  because those fields are in the manifest. The release guard refuses an existing
  tag whose assets differ, so such a rebuild needs a new revision before it can be
  published.

## Making a release

1. Merge the recipe to `dev` (the release guard refuses a `recipe_commit` that
   is not on `origin/dev`).
2. Build from a clean checkout of the `origin/dev` commit you will publish, into
   an empty `dist/<id>`. The build records that commit itself.
3. Dry run: `python -c "from pathlib import Path; from feeds.release import publish; publish(Path('dist/<id>'))"`
   (it prints the `gh release create` command it would run). A dry run is not
   offline: it runs `git fetch origin dev` in this checkout (after checking that
   `origin` is `chd-atlas/feeds`), `gh release view` and, when the tag exists,
   `gh release download`. It only withholds `gh release create`.
4. The owner publishes with `dry_run=False`. A tag that exists with identical
   assets is `unchanged`; one with different assets is refused. The release notes
   state the source's licence and redistribution terms, taken from the manifest.

## Licence

The `LICENSE` in this repository (Apache-2.0) covers the code rather than the
data of the packages it publishes. Each package's data carries its source's
terms, stated in its `manifest.json` (`source.licence`, `source.redistribution`,
`licence_evidence`) and repeated on its release page. For example, HGNC is CC0 and
Sierant et al. 2025 is CC BY 4.0; Farah et al. 2024's data state no licence, and
are treated as attributable to the article (see that source's README).
