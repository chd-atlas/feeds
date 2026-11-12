# gencc

GenCC submissions export (CC0-1.0) as the `gencc` package: table `gencc_submissions`
v1, `ids: hgnc_id`, `coverage: every_gene`. GenCC states no release, so the version
is the input's sha256[:8] and `source.version` is `retrieved <retrieved_on>`.

Every submission is kept; the key is `(gene, disease, moi, submitter, sgc_id)` (see
`recipe.py` for why `sgc_id` is in it). The export is quoted TSV, parsed with polars,
never by splitting lines.

Two entry points in `recipe.py`:

- `build(tsv, retrieved_on=..., revision=..., built_by=..., out=...)` parses a live
  download from `https://search.thegencc.org/download/action/submissions-export-tsv?format=new`.
- `bootstrap(mirror, atlas_commit=..., ...)` packages the atlas's committed
  `mirrors/gencc_submissions.tsv` once (decision F21), because the 2026-08-03 export
  cannot be fetched again. It accepts only the file with sha256 `88c77272...` (atlas
  4b1f8f6, 30,410 rows, 6,113 genes); the package version is `88c77272`.

## Build

Live (`?format=new` is mandatory; the quota is 20 requests per IP per day; GenCC states
no release, so the sha256 is whatever you got):

```bash
mkdir -p local/gencc
curl -L -o local/gencc/gencc-submissions.tsv "https://search.thegencc.org/download/action/submissions-export-tsv?format=new"
```

```python
from datetime import date
from pathlib import Path
from sources.gencc import recipe

recipe.build(
    Path("local/gencc/gencc-submissions.tsv"),
    retrieved_on=date.today(),
    revision=1,
    built_by={"host": "<machine>"},
    out=Path("dist/gencc"),
)
```

Bootstrap, from the atlas mirror at the pinned commit (`recipe.ATLAS_COMMIT`,
`4b1f8f68cbc30beaf932fbbc58e3617ebe401b6d`; any other `atlas_commit` is refused, because
the manifest's input URL is a permanent permalink to it):

```bash
mkdir -p local/gencc
git -C <atlas> show 4b1f8f6:mirrors/gencc_submissions.tsv > local/gencc/gencc_submissions.tsv
```

```python
recipe.bootstrap(
    Path("local/gencc/gencc_submissions.tsv"),
    atlas_commit=recipe.ATLAS_COMMIT,
    revision=1,
    built_by={"host": "<machine>"},
    out=Path("dist/gencc"),
)
```

## Measurements

- Live parity, 2026-11-12: `read_rows` over the 2026-11-11 export gives 30,376
  records (46,968 physical lines) and `dropped == []`; rendered as TSV it is
  `cmp`-equal to the atlas's `scripts/convert_gencc.py` (archive of atlas 4b1f8f6)
  on the same file.
- A bootstrap build of the real mirror into a scratch directory succeeded; its table
  is byte-identical to the mirror.

## Tests

`pytest sources/gencc`: synthetic export and mirror, no network.
