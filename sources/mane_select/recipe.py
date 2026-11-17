"""MANE v1.5 Select transcripts and their exons (+-50 bp of GRCh38) as the ``mane_select`` package.

Two tables, whose schemas the atlas exports: ``mane_transcripts`` (one row per gene: the
transcript and protein sequences and the coding bounds) and ``mane_exons`` (one row per
exon: its transcript and genomic coordinates and the forward-strand genome window that
holds the exon and ``FLANK`` bases either side). The atlas checks the two against each
other (MAN001, MAN002) and against its gene registry (MAN003, MAN004); this recipe runs
MAN001 and MAN002 itself before writing, so a transcript that would fail them is listed
in ``dropped`` with its reason instead of reaching the atlas.

INPUTS
------
The summary, the GTF, the RNA and protein FASTA of MANE v1.5, and the GRCh38.p14 genome
from NCBI. All five are pinned by sha256 (``PINNED_SHA256``); any other bytes are refused
before a file is written.

CODING BOUNDS
-------------
``cds_start`` is the transcript position of the first ``CDS`` base and ``cds_end`` that of
the last ``CDS`` base plus 3, because the GTF's ``CDS`` feature excludes the stop codon.
``start_codon`` and ``stop_codon`` features are not used: they are not one per transcript
(a codon split by an intron, or absent), whereas every coding transcript has CDS. MAN002
then confirms the three bases after the CDS really are a stop.

WHAT IS DROPPED, AND WHY
------------------------
* A MANE Select transcript that is non-coding (an ``NR_`` accession): it has no protein.
* A MANE Select transcript of a gene with no HGNC id (``LOC...`` genes): the package is
  keyed by HGNC id, and the atlas's gene registry is.
* A MANE Select transcript the summary places on a non-primary accession (``NT_``,
  ``NW_``): the atlas's chromosome vocabulary is 1-22, X, Y, MT, and a position on an
  unlocalised scaffold would be a guess.
* A transcript whose exon cores do not reproduce its RNA sequence (MAN001) or whose
  coding sequence does not translate to its protein (MAN002). MAN002 allows ``U``
  against ``TGA`` (selenocysteine) and ``M`` at residue 1 against a non-ATG start,
  exactly as the atlas does; anything else is dropped, never repaired.
* A transcript whose sequence or identifiers do not fit the table's patterns.
More than ``MAX_FAILURES`` such failures raises ``PackageRefused``: that many would mean
the parser is wrong, not MANE.

AMBIGUITY CODES
---------------
``window_sequence`` allows ``ACGTN``. An IUPAC code in a FLANK is written as ``N`` (and
counted, in the README); one inside an exon core cannot match the transcript, so that
transcript is dropped. Windows near a chromosome end are clipped; ``window_start`` and
``window_end`` are the truth.
"""

from __future__ import annotations

import gzip
import hashlib
import re
from collections import defaultdict
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any, Final

from feeds.package import PackageRefused, TableOut, write_package
from feeds.provenance import CommitOf, current_commit, stamp
from feeds.source import load_source
from feeds.version import package_version

HERE = Path(__file__).resolve().parent
SOURCE: Final = "mane_select"
UPSTREAM_RELEASE: Final = "1.5"
FLANK: Final = 50
MAX_FAILURES: Final = 50
GENOME_NAME: Final = "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
_BASE_URL: Final = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/"
SUMMARY: Final = "MANE.GRCh38.v1.5.summary.txt.gz"
GTF: Final = "MANE.GRCh38.v1.5.refseq_genomic.gtf.gz"
RNA: Final = "MANE.GRCh38.v1.5.refseq_rna.fna.gz"
PROTEIN: Final = "MANE.GRCh38.v1.5.refseq_protein.faa.gz"
DOWNLOAD_URLS: Final[dict[str, str]] = {
    SUMMARY: _BASE_URL + SUMMARY,
    GTF: _BASE_URL + GTF,
    RNA: _BASE_URL + RNA,
    PROTEIN: _BASE_URL + PROTEIN,
    GENOME_NAME: (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
        "GCF_000001405.40_GRCh38.p14/" + GENOME_NAME
    ),
}
PINNED_SHA256: dict[str, str] = {
    SUMMARY: "d10ace2720681a3b2e0eefd9da4f551274a6b4141ac9bfd6a2565dfb6e9ad55c",
    GTF: "59e5ff1c794a70de777e553ab67062b316578d41287cd898ee1da2827527969a",
    RNA: "7280512e57fd8599ab4e3a4a2e809b7f1780c3b9e30a8ca78ca937a7d8dee98c",
    PROTEIN: "b12307f98481b6d8c4e17d9ecd1856ee1f64a63db6507ff7529f669e26f6f372",
    GENOME_NAME: "11912a45a545bf01a10b2a7f10eb7a42924436b4d19b476b1899834fb7ba74a3",
}
# From the release's README_versions.txt, read 2026-11-17.
FILE_VERSION: Final = (
    "MANE 1.5; NCBI RefSeq Annotation Release GCF_000001405.40-RS_2025_08; Ensembl 116"
)

TRANSCRIPT_COLUMNS: Final = (
    "hgnc_id",
    "refseq_transcript",
    "refseq_protein",
    "ensembl_transcript",
    "ensembl_protein",
    "chrom",
    "strand",
    "cds_start",
    "cds_end",
    "transcript_sequence",
    "protein_sequence",
    "source",
)
EXON_COLUMNS: Final = (
    "hgnc_id",
    "refseq_transcript",
    "exon",
    "transcript_start",
    "transcript_end",
    "genomic_start",
    "genomic_end",
    "window_start",
    "window_end",
    "window_sequence",
    "source",
)
# The atlas's column patterns (tables.py MANE_TRANSCRIPTS / MANE_EXONS), so a value that
# would refuse the whole package drops its one transcript instead.
_PATTERNS: Final = {
    "hgnc_id": re.compile(r"HGNC:\d+"),
    "refseq_transcript": re.compile(r"NM_\d+\.\d+"),
    "refseq_protein": re.compile(r"NP_\d+\.\d+"),
    "ensembl_transcript": re.compile(r"ENST\d{11}\.\d+"),
    "ensembl_protein": re.compile(r"ENSP\d{11}\.\d+"),
    "transcript_sequence": re.compile(r"[ACGT]+"),
    "protein_sequence": re.compile(r"[ACDEFGHIKLMNPQRSTVWYU]+"),
}

_PRIMARY = {f"NC_{n:06d}": str(n) for n in range(1, 23)}
_PRIMARY["NC_000023"] = "X"
_PRIMARY["NC_000024"] = "Y"
_PRIMARY["NC_012920"] = "MT"

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE: Final[dict[str, str]] = {
    a + b + c: _AMINO[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
_NOT_BASE = re.compile(r"[^ACGTN]")
_EXAMPLES = 5
_FAULT_REASON: Final = {
    "MAN001": "MAN001: the exon cores do not reproduce the RNA sequence",
    "MAN002": "MAN002: the coding sequence does not translate to the protein",
}


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class _Tx:
    hgnc_id: str
    nuc: str
    prot: str
    ens_nuc: str
    ens_prot: str
    nc: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)  # genomic, transcript order
    cds: list[tuple[int, int]] = field(default_factory=list)  # genomic
    cds_start: int = 0
    cds_end: int = 0
    spans: list[tuple[int, int]] = field(default_factory=list)  # transcript start/end per exon
    windows: list[tuple[int, int, str]] = field(default_factory=list)


class _Dropped:
    """Dropped transcripts grouped by reason, with up to five examples each."""

    def __init__(self) -> None:
        self.by_reason: dict[str, list[str]] = defaultdict(list)
        self.failures = (
            0  # every drop except the by-design ones (non-coding, no HGNC id, non-primary)
        )

    def add(self, reason: str, example: str, *, failure: bool = True) -> None:
        self.by_reason[reason].append(example)
        self.failures += failure

    def entries(self) -> list[dict[str, Any]]:
        return [
            {"reason": reason, "count": len(items), "examples": sorted(items)[:_EXAMPLES]}
            for reason, items in sorted(self.by_reason.items())
        ]


def _lines(path: Path) -> Iterator[str]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            yield line.rstrip("\r\n")


def _fasta(path: Path) -> dict[str, str]:
    """``{accession: sequence}``, the accession being the first token of the header."""
    records: dict[str, list[str]] = {}
    current: list[str] | None = None
    for line in _lines(path):
        if line.startswith(">"):
            current = records.setdefault(line[1:].split(None, 1)[0], [])
        elif current is not None:
            current.append(line.strip())
    return {name: "".join(parts).upper() for name, parts in records.items()}


def _summary(path: Path, dropped: _Dropped) -> tuple[list[_Tx], int]:
    """``(MANE Select transcripts on a primary accession, MANE Select rows read)``."""
    lines = _lines(path)
    header = next(lines, "").lstrip("#").split("\t")
    want = (
        "HGNC_ID",
        "RefSeq_nuc",
        "RefSeq_prot",
        "Ensembl_nuc",
        "Ensembl_prot",
        "MANE_status",
        "GRCh38_chr",
        "chr_strand",
    )
    missing = [c for c in want if c not in header]
    if missing:
        raise PackageRefused(f"{path.name}: no column {', '.join(missing)}")
    at = {c: header.index(c) for c in want}
    txs: list[_Tx] = []
    read = 0
    seen: set[str] = set()
    for line in lines:
        cells = line.split("\t")
        if not line or len(cells) <= max(at.values()) or cells[at["MANE_status"]] != "MANE Select":
            continue
        read += 1
        hgnc = cells[at["HGNC_ID"]]
        nc = cells[at["GRCh38_chr"]]
        nuc = cells[at["RefSeq_nuc"]]
        if not nuc.startswith("NM_"):
            dropped.add(
                "MANE Select non-coding transcript (NR_), which has no protein to frame",
                f"{hgnc or 'no HGNC id'} {nuc}",
                failure=False,
            )
            continue
        if not _PATTERNS["hgnc_id"].fullmatch(hgnc):
            dropped.add(
                "MANE Select transcript of a gene with no HGNC id (the table is keyed by HGNC id)",
                f"{nuc} ({cells[at['HGNC_ID']] or 'blank'})",
                failure=False,
            )
            continue
        if hgnc in seen:
            raise PackageRefused(f"{path.name}: {hgnc} has more than one MANE Select row")
        seen.add(hgnc)
        chrom = _PRIMARY.get(nc.split(".")[0])
        if chrom is None:
            dropped.add(
                "MANE Select transcript placed on a non-primary accession "
                "(the atlas's chromosome vocabulary is 1-22, X, Y, MT)",
                f"{hgnc} {nuc} on {nc}",
                failure=False,
            )
            continue
        txs.append(
            _Tx(
                hgnc,
                nuc,
                cells[at["RefSeq_prot"]],
                cells[at["Ensembl_nuc"]],
                cells[at["Ensembl_prot"]],
                nc,
                chrom,
                cells[at["chr_strand"]],
            )
        )
    return txs, read


_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]*)"')
_SELECT_TAG = 'tag "MANE Select";'


def _read_gtf(path: Path, txs: dict[str, _Tx]) -> None:
    """Attach the exons and CDS of each Select transcript, on the chromosome the summary names."""
    for line in _lines(path):
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9 or cols[2] not in ("exon", "CDS") or _SELECT_TAG not in cols[8]:
            continue
        match = _TRANSCRIPT_ID.search(cols[8])
        tx = txs.get(match.group(1)) if match else None
        if tx is None or cols[0] != f"chr{tx.chrom}" or cols[6] != tx.strand:
            continue
        span = (int(cols[3]), int(cols[4]))
        (tx.exons if cols[2] == "exon" else tx.cds).append(span)


def _lay_out(tx: _Tx) -> str | None:
    """Order the exons, cumulate their transcript spans and place the CDS. None when it fits."""
    if not tx.exons:
        return "no exon in the GTF on the summary's chromosome and strand"
    if not tx.cds:
        return "no CDS in the GTF"
    minus = tx.strand == "-"
    tx.exons.sort(reverse=minus)
    at = 1
    for start, end in tx.exons:
        tx.spans.append((at, at + end - start))
        at += end - start + 1

    def locate(pos: int) -> int | None:
        for (start, end), (t_start, _) in zip(tx.exons, tx.spans, strict=True):
            if start <= pos <= end:
                return t_start + (end - pos if minus else pos - start)
        return None

    placed = [locate(p) for span in tx.cds for p in span]
    if None in placed:
        return "a CDS base lies outside the exons"
    bases = [p for p in placed if p is not None]
    tx.cds_start, tx.cds_end = min(bases), max(bases) + 3  # the GTF's CDS excludes the stop
    return None


def _cut_windows(txs: list[_Tx], genome: Path, dropped: _Dropped) -> list[_Tx]:
    """Cut each exon's window from the genome, streaming it once, one chromosome at a time."""
    need: dict[str, list[_Tx]] = defaultdict(list)
    for tx in txs:
        need[tx.nc].append(tx)
    done: set[str] = set()

    def cut(accession: str, parts: list[bytes]) -> None:
        sequence = b"".join(parts).decode("ascii").upper()
        length = len(sequence)
        for tx in need[accession]:
            for start, end in tx.exons:
                w_start, w_end = max(1, start - FLANK), min(length, end + FLANK)
                tx.windows.append((w_start, w_end, sequence[w_start - 1 : w_end]))
        done.add(accession)

    current: str | None = None
    parts: list[bytes] = []
    with gzip.open(genome, "rb") as handle:
        for raw in handle:
            if raw.startswith(b">"):
                if current is not None:
                    cut(current, parts)
                token = raw[1:].split(None, 1)[0].decode("ascii")
                current = token if token in need and token not in done else None
                parts = []
            elif current is not None:
                parts.append(raw.strip())
    if current is not None:
        cut(current, parts)
    kept: list[_Tx] = []
    for tx in txs:
        if tx.nc not in done:
            dropped.add(f"chromosome {tx.nc} is not in the genome FASTA", f"{tx.hgnc_id} {tx.nuc}")
        else:
            kept.append(tx)
    return kept


def _fault(tx: _Tx, sequence: str, protein: str) -> tuple[str, str] | None:
    """``(code, reason)`` for the first way the transcript fails MAN001 or MAN002, else None.

    Written to agree with the atlas's ``validate/mane.py`` (MAN001 ``_frame_fault``, MAN002
    ``_cds_fault``), including MAN002's two allowances: ``U`` against ``TGA`` and ``M`` at
    residue 1 against a codon other than ``ATG``.
    """
    pieces: list[str] = []
    for (g_start, g_end), (w_start, _, window) in zip(tx.exons, tx.windows, strict=True):
        core = window[g_start - w_start : g_end - w_start + 1]
        if _NOT_BASE.search(core) or "N" in core:
            return "MAN001", "an exon core holds a base other than A, C, G or T"
        pieces.append(_revcomp(core) if tx.strand == "-" else core)
    joined = "".join(pieces)
    if joined != sequence:
        first = next(
            (i + 1 for i, (a, b) in enumerate(zip(joined, sequence, strict=False)) if a != b),
            min(len(joined), len(sequence)) + 1,
        )
        return "MAN001", f"exon sequences differ from the RNA sequence at position {first}"
    if not 1 <= tx.cds_start <= tx.cds_end <= len(sequence):
        return "MAN002", f"cds {tx.cds_start}-{tx.cds_end} is not inside the sequence"
    cds = sequence[tx.cds_start - 1 : tx.cds_end]
    if len(cds) % 3:
        return "MAN002", f"cds length {len(cds)} is not a multiple of 3"
    codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
    translated = "".join(CODON_TABLE.get(c, "?") for c in codons)
    if translated[-1] != "*":
        return "MAN002", "cds does not end in a stop codon"
    body = translated[:-1]
    if len(body) != len(protein):
        return "MAN002", f"cds encodes {len(body)} residues but the protein has {len(protein)}"
    for index, (got, want, codon) in enumerate(zip(body, protein, codons[:-1], strict=True)):
        if got == want:
            continue
        if want == "U" and codon == "TGA":
            continue
        if index == 0 and want == "M" and codon != "ATG":
            continue
        return "MAN002", f"residue {index + 1} is {want} in the protein but {codon} encodes {got}"
    return None


def _pattern_fault(tx: _Tx, sequence: str, protein: str) -> str | None:
    values = {
        "hgnc_id": tx.hgnc_id,
        "refseq_transcript": tx.nuc,
        "refseq_protein": tx.prot,
        "ensembl_transcript": tx.ens_nuc,
        "ensembl_protein": tx.ens_prot,
        "transcript_sequence": sequence,
        "protein_sequence": protein,
    }
    for name, value in values.items():
        if not _PATTERNS[name].fullmatch(value):
            return name
    return None


@dataclass(frozen=True)
class Built:
    transcripts: list[tuple[str, ...]]
    exons: list[tuple[str, ...]]
    dropped: list[dict[str, Any]]
    rows_read: int
    flank_n: int  # flank bases already N in the genome
    flank_iupac: int  # flank bases with another IUPAC code, written as N


def read_rows(summary: Path, gtf: Path, rna: Path, protein: Path, genome: Path) -> Built:
    """The package rows and ``dropped`` from the five files. Pure: no sha256 pin is checked."""
    dropped = _Dropped()
    txs, read = _summary(summary, dropped)
    by_nuc = {tx.nuc: tx for tx in txs}
    _read_gtf(gtf, by_nuc)
    rna_by = _fasta(rna)
    protein_by = _fasta(protein)
    laid: list[_Tx] = []
    for tx in txs:
        why = _lay_out(tx)
        if why is not None:
            dropped.add(why, f"{tx.hgnc_id} {tx.nuc}")
        else:
            laid.append(tx)
    cut = _cut_windows(laid, genome, dropped)
    rows_t: list[tuple[str, ...]] = []
    rows_e: list[tuple[str, ...]] = []
    flank_n = flank_iupac = 0
    for tx in cut:
        sequence, prot = rna_by.get(tx.nuc), protein_by.get(tx.prot)
        if sequence is None or prot is None:
            which = "RNA" if sequence is None else "protein"
            dropped.add(f"no {which} sequence in the MANE FASTA", f"{tx.hgnc_id} {tx.nuc}")
            continue
        bad = _pattern_fault(tx, sequence, prot)
        if bad is not None:
            dropped.add(f"{bad} does not fit the table's pattern", f"{tx.hgnc_id} {tx.nuc}")
            continue
        fault = _fault(tx, sequence, prot)
        if fault is not None:
            code, reason = fault
            dropped.add(_FAULT_REASON[code], f"{tx.hgnc_id} {tx.nuc}: {reason}")
            continue
        tx_exons: list[tuple[str, ...]] = []
        for number, ((g_start, g_end), (t_start, t_end), (w_start, w_end, window)) in enumerate(
            zip(tx.exons, tx.spans, tx.windows, strict=True), start=1
        ):
            lead = g_start - w_start
            tail = g_end - w_start + 1
            flanks = window[:lead], window[tail:]
            flank_n += sum(f.count("N") for f in flanks)
            flank_iupac += sum(len(_NOT_BASE.findall(f)) for f in flanks)
            window = _NOT_BASE.sub("N", window[:lead]) + window[lead:tail]
            window += _NOT_BASE.sub("N", flanks[1])
            tx_exons.append(
                (
                    tx.hgnc_id,
                    tx.nuc,
                    str(number),
                    str(t_start),
                    str(t_end),
                    str(g_start),
                    str(g_end),
                    str(w_start),
                    str(w_end),
                    window,
                    SOURCE,
                )
            )
        rows_e.extend(tx_exons)
        rows_t.append(
            (
                tx.hgnc_id,
                tx.nuc,
                tx.prot,
                tx.ens_nuc,
                tx.ens_prot,
                tx.chrom,
                tx.strand,
                str(tx.cds_start),
                str(tx.cds_end),
                sequence,
                prot,
                SOURCE,
            )
        )
    if dropped.failures > MAX_FAILURES:
        raise PackageRefused(
            f"{dropped.failures} MANE Select transcripts failed (more than {MAX_FAILURES}): "
            "that many means the parser is wrong, not MANE. First reasons: "
            + "; ".join(f"{r} x{len(i)}" for r, i in sorted(dropped.by_reason.items())[:3])
        )
    rows_t.sort(key=lambda r: r[0])
    rows_e.sort(key=lambda r: (r[0], int(r[2])))
    return Built(rows_t, rows_e, dropped.entries(), read, flank_n, flank_iupac)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def build(
    summary: Path,
    gtf: Path,
    rna: Path,
    protein: Path,
    genome: Path,
    *,
    retrieved_on: date,
    revision: int,
    built_by: dict[str, Any],
    out: Path,
    commit_of: CommitOf = current_commit,
) -> Path:
    """Write the ``mane_select`` package from the five pinned downloads into the empty ``out``.

    ``revision`` has no default (F13): 1 for the first build of these inputs, 2, 3, ... for
    a recipe change or a rebuild whose ``manifest.json`` bytes differ.
    """
    built_by, _ = stamp(built_by, commit_of)
    files = {SUMMARY: summary, GTF: gtf, RNA: rna, PROTEIN: protein, GENOME_NAME: genome}
    digests = {name: _sha256(path) for name, path in files.items()}
    for name, sha in digests.items():
        if sha != PINNED_SHA256[name]:
            raise PackageRefused(
                f"{name} has sha256 {sha}, not the pinned {PINNED_SHA256[name]}; "
                "this recipe packages MANE v1.5 and GRCh38.p14 only"
            )
    built = read_rows(summary, gtf, rna, protein, genome)
    source = load_source(HERE / "source.yaml")
    stamped = retrieved_on.isoformat()
    transcripts: list[Sequence[str | None]] = list(built.transcripts)
    exons: list[Sequence[str | None]] = list(built.exons)
    return write_package(
        out,
        package=SOURCE,
        version=package_version(
            upstream_release=UPSTREAM_RELEASE, input_sha256=digests[SUMMARY], revision=revision
        ),
        ids="hgnc_id",
        coverage="every_gene",
        tables=[
            TableOut("mane_transcripts", 1, TRANSCRIPT_COLUMNS, transcripts, "hgnc_id"),
            TableOut("mane_exons", 1, EXON_COLUMNS, exons, "hgnc_id"),
        ],
        inputs=[
            {
                "name": name,
                "url": DOWNLOAD_URLS[name],
                "sha256": digests[name],
                "retrieved_on": stamped,
                "upstream_release": (
                    "GCF_000001405.40 (GRCh38.p14)" if name == GENOME_NAME else UPSTREAM_RELEASE
                ),
            }
            for name in (SUMMARY, GTF, RNA, PROTEIN, GENOME_NAME)
        ],
        input_counts={"rows": built.rows_read},
        dropped=built.dropped,
        built_by=built_by,
        source={
            "name": source["name"],
            "version": f"MANE v{UPSTREAM_RELEASE}",
            "file_version": FILE_VERSION,
            "retrieved_on": stamped,
            "url": source["url"],
            "licence": source["licence"],
            "redistribution": source["redistribution"],
        },
        licence_evidence=source["licence_evidence"],
    )
