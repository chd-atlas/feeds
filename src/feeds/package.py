"""Write a source package: gzipped TSV tables plus ``manifest.json``.

Every refusal happens before any byte reaches ``out_dir``. The files are
byte-reproducible: sorted-key manifest, gzip with ``mtime=0`` and no file name.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import shutil
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import jsonschema

from feeds.conform import check_table
from feeds.version import SEGMENT

VENDORED = Path(__file__).resolve().parents[2] / "schemas" / "atlas"


# The source that owns each flat table without a source column (atlas ING014): another
# source's package would replace all of it. The atlas's second ING014 branch, a shard
# another committed manifest carries, depends on its repository and is not mirrored.
FLAT_TABLE_OWNERS = {"clingen_validity": "clingen", "gencc_submissions": "gencc"}
HGNC = "hgnc"


class PackageRefused(Exception):
    """The package would be misread or refused by the atlas; nothing was written."""


@dataclass(frozen=True)
class TableOut:
    name: str
    format_version: int
    columns: tuple[str, ...]
    rows: Sequence[Sequence[str | None]]
    gene_column: str | None
    shard: str | None = None


def _encode(table: TableOut) -> bytes:
    lines = ["\t".join(table.columns)]
    for number, row in enumerate(table.rows, start=1):
        for at, cell in enumerate(row):
            if cell is None:
                continue
            column = table.columns[at] if at < len(table.columns) else f"#{at + 1}"
            if "\t" in cell or "\n" in cell or "\r" in cell:
                raise PackageRefused(
                    f"table '{table.name}' row {number}, column '{column}': a cell holds a tab, "
                    "newline or carriage return"
                )
            if '"' in cell:
                raise PackageRefused(
                    f"table '{table.name}' row {number}, column '{column}': a cell holds a "
                    "double quote (the atlas refuses a leading one as ING015 and misreads "
                    "an inner one as ING011)"
                )
        lines.append("\t".join("" if cell is None else cell for cell in row))
    try:
        raw = ("\n".join(lines) + "\n").encode("utf-8")
    except UnicodeEncodeError as exc:
        raise PackageRefused(
            f"table '{table.name}' holds text that is not valid UTF-8: {exc}"
        ) from exc
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer, mtime=0) as handle:
        handle.write(raw)
    return buffer.getvalue()


def _load(path: Path) -> dict[str, Any]:
    loaded: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return loaded


def write_package(
    out_dir: Path,
    *,
    package: str,
    version: str,
    ids: str,
    coverage: str | dict[str, Any],
    tables: Sequence[TableOut],
    inputs: Sequence[dict[str, Any]],
    input_counts: dict[str, Any],
    dropped: Sequence[dict[str, Any]],
    built_by: dict[str, Any],
    source: dict[str, Any],
    licence_evidence: str,
    schema_dir: Path = VENDORED,
) -> Path:
    if not re.fullmatch(SEGMENT, version):
        raise PackageRefused(f"version '{version}' is not a valid path segment")
    if out_dir.exists() and (not out_dir.is_dir() or any(out_dir.iterdir())):
        raise PackageRefused(f"{out_dir} exists and is not empty; never mix two packages")
    if not licence_evidence.strip():
        raise PackageRefused("licence_evidence must not be blank")
    names = [t.name for t in tables]
    # The atlas's ING009: `genes` is the registry, and only the hgnc package ships it.
    if package != HGNC and "genes" in names:
        raise PackageRefused(
            f"package '{package}' ships the genes table; only '{HGNC}' may (atlas ING009)"
        )
    if package == HGNC and "genes" not in names:
        raise PackageRefused("the hgnc package must ship the genes table (atlas ING009)")
    if package == HGNC and ids != "hgnc_id":
        raise PackageRefused(f"the hgnc package's ids must be hgnc_id, not '{ids}' (atlas ING009)")
    for name in names:
        owner = FLAT_TABLE_OWNERS.get(name)
        if owner is not None and owner != package:
            raise PackageRefused(
                f"table '{name}' belongs to source '{owner}'; package '{package}' would "
                "replace all of it (atlas ING014)"
            )
    repeated = sorted({n for n in names if names.count(n) > 1})
    if repeated:
        raise PackageRefused(f"table names repeated: {', '.join(repeated)}")
    entries: list[dict[str, Any]] = []
    blobs: dict[str, bytes] = {}
    for table in tables:
        if table.shard is not None and not re.fullmatch(SEGMENT, table.shard):
            raise PackageRefused(f"shard '{table.shard}' is not a valid path segment")
        schema_path = schema_dir / "tables" / f"{table.name}.v{table.format_version}.json"
        if not schema_path.is_file():
            raise PackageRefused(f"no vendored schema for {table.name} v{table.format_version}")
        schema = _load(schema_path)
        if table.gene_column != schema.get("geneColumn"):
            raise PackageRefused(
                f"table '{table.name}': gene_column {table.gene_column!r} differs from the "
                f"schema's {schema.get('geneColumn')!r}"
            )
        layout = schema.get("layout", {})
        if "sharded" in layout and table.shard is None:
            raise PackageRefused(f"table '{table.name}' is sharded and gives no shard")
        if "flat" in layout and table.shard is not None:
            raise PackageRefused(f"table '{table.name}' is flat and gives a shard")
        source_column = schema.get("sourceColumn")
        if source_column is not None:
            at_source = table.columns.index(source_column) if source_column in table.columns else -1
            foreign = sorted(
                {
                    repr(row[at_source])
                    for row in table.rows
                    if at_source >= 0 and row[at_source] != package
                }
            )
            if foreign:
                raise PackageRefused(
                    f"table '{table.name}' has rows whose {source_column} is not '{package}': "
                    f"{foreign[:5]}"
                )
        errors = check_table(table.columns, table.rows, schema, ids=ids)
        if errors:
            shown = "; ".join(errors[:5])
            raise PackageRefused(f"table '{table.name}' does not conform ({len(errors)}): {shown}")
        blob = _encode(table)
        file = f"{table.name}.{table.shard}.tsv.gz" if table.shard else f"{table.name}.tsv.gz"
        blobs[file] = blob
        genes: int | None = None
        if table.gene_column is not None:
            at = table.columns.index(table.gene_column)
            stripped = {(row[at] or "").strip() for row in table.rows}
            genes = len(stripped - {""})
        entry: dict[str, Any] = {
            "name": table.name,
            "format_version": table.format_version,
            "file": file,
            "sha256": hashlib.sha256(blob).hexdigest(),
            "rows": len(table.rows),
            "genes": genes,
        }
        if table.shard is not None:
            entry["shard"] = table.shard
        entries.append(entry)
    entries.sort(key=lambda e: (e["name"], e.get("shard") or ""))
    manifest: dict[str, Any] = {
        "package": package,
        "version": version,
        "ids": ids,
        "coverage": coverage,
        "tables": entries,
        "inputs": list(inputs),
        "input_counts": input_counts,
        "dropped": list(dropped),
        "built_by": built_by,
        "source": source,
        "licence_evidence": licence_evidence,
    }
    manifest_schema = _load(schema_dir / "package_manifest.schema.json")
    validator_class = jsonschema.validators.validator_for(
        manifest_schema, default=jsonschema.Draft202012Validator
    )
    problems = sorted(
        validator_class(manifest_schema, format_checker=validator_class.FORMAT_CHECKER).iter_errors(
            manifest
        ),
        key=lambda e: list(e.path),
    )
    if problems:
        shown = "; ".join(f"{'/'.join(map(str, e.path))}: {e.message}" for e in problems[:5])
        raise PackageRefused(f"manifest does not satisfy the atlas's schema: {shown}")
    try:
        text = (json.dumps(manifest, sort_keys=True, indent=2, ensure_ascii=False) + "\n").encode()
    except UnicodeEncodeError as exc:
        raise PackageRefused(f"manifest holds text that is not valid UTF-8: {exc}") from exc
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{out_dir.name}.", dir=out_dir.parent))
    try:
        for file, blob in blobs.items():
            (staging / file).write_bytes(blob)
        (staging / "manifest.json").write_bytes(text)
        staging.replace(out_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return out_dir
