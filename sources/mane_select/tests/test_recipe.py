"""The MANE Select recipe on a synthetic genome: no network, no committed MANE file.

``World`` writes a made-up genome and the four MANE files for the genes it is given. Each
transcript is written first (5' UTR, a coding sequence back-translated from a protein, a
stop codon, 3' UTR) and then laid into a chromosome at the given positions, on the given
strand, so exon cores, the RNA FASTA and the GTF agree by construction and every
disagreement a test wants is introduced on purpose. Codons are written out below, not
inverted from the code under test.
"""

import gzip
import hashlib
import json
import os
import random
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

import pytest

from feeds.package import PackageRefused
from sources.mane_select import recipe

BUILT_BY = {"host": "login01"}
RECIPE_COMMIT = "a" * 40
CODON = {
    "M": "ATG",
    "K": "AAA",
    "T": "ACC",
    "A": "GCC",
    "R": "CGG",
    "S": "TCC",
    "L": "CTG",
    "G": "GGC",
    "W": "TGG",
    "C": "TGC",
    "U": "TGA",
}
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


@dataclass
class Gene:
    hgnc: str
    nuc: str
    protein: str  # what the protein FASTA says
    nc: str = "NC_000012.12"
    strand: str = "+"
    utr5: str = "GGAATTCCGG"
    utr3: str = "CCTTAAGGCCTT"
    cuts: tuple[int, ...] = ()  # transcript lengths of every exon but the last
    starts: tuple[int, ...] = ()  # genomic start of each exon, in genome order
    coding: str | None = None  # codons written into the transcript, when not the protein's
    status: str = "MANE Select"
    extra_gtf: list[str] = field(default_factory=list)

    @property
    def cds(self) -> str:
        codons = self.coding or "".join(CODON[a] for a in self.protein)
        return codons + "TAA"

    @property
    def transcript(self) -> str:
        return self.utr5 + self.cds + self.utr3

    def exon_lengths(self) -> list[int]:
        total, lengths, last = len(self.transcript), [], 0
        for cut in (*self.cuts, total):
            lengths.append(cut - last)
            last = cut
        return lengths


class World:
    """A genome with chromosome ``NC_000012.12`` (and optionally others) and the MANE files."""

    def __init__(self, tmp: Path, *, length: int = 3000, seed: int = 5) -> None:
        self.tmp, self.length = tmp, length
        rng = random.Random(seed)
        self.chrom: dict[str, list[str]] = {}
        self.rng = rng
        self.genes: list[Gene] = []
        self.add_chrom("NC_000012.12")

    def add_chrom(self, nc: str, length: int | None = None) -> None:
        # Soft-masked, as the real genome is: the recipe must upper-case it.
        self.chrom[nc] = [self.rng.choice("acgtACGT") for _ in range(length or self.length)]

    def add(self, gene: Gene) -> Gene:
        self.genes.append(gene)
        if gene.nc not in self.chrom:
            self.add_chrom(gene.nc)
        if not gene.starts:
            return gene
        seq, at = gene.transcript, 0
        pieces = []
        for length in gene.exon_lengths():
            pieces.append(seq[at : at + length])
            at += length
        ordered = pieces[::-1] if gene.strand == "-" else pieces
        for start, piece in zip(gene.starts, ordered, strict=True):
            core = revcomp(piece) if gene.strand == "-" else piece
            self.chrom[gene.nc][start - 1 : start - 1 + len(core)] = list(core.lower())
        return gene

    def spans(self, gene: Gene) -> list[tuple[int, int]]:
        """Genomic (start, end) of each exon in transcript order."""
        lengths = gene.exon_lengths()
        by_genome = lengths[::-1] if gene.strand == "-" else lengths
        spans = [(s, s + n - 1) for s, n in zip(gene.starts, by_genome, strict=True)]
        return spans[::-1] if gene.strand == "-" else spans

    def write(self) -> dict[str, Path]:
        d = self.tmp / "in"
        d.mkdir(exist_ok=True)
        summary = ["#NCBI_GeneID\tEnsembl_Gene\tHGNC_ID\tsymbol\tname\tRefSeq_nuc\tRefSeq_prot"]
        summary[0] += "\tEnsembl_nuc\tEnsembl_prot\tMANE_status\tGRCh38_chr\tchr_start\tchr_end"
        summary[0] += "\tchr_strand"
        gtf = ["#gtf-version 2.2"]
        rna, prot = [], []
        for n, g in enumerate(self.genes, start=1):
            num = f"{n:011d}"
            summary.append(
                f"GeneID:{n}\tENSG{num}.1\t{g.hgnc}\tG{n}\tname\t{g.nuc}"
                f"\t{g.nuc.replace('NM', 'NP')}"
                f"\tENST{num}.1\tENSP{num}.1\t{g.status}\t{g.nc}\t1\t2\t{g.strand}"
            )
            rna += [f">{g.nuc} Homo sapiens test, mRNA", g.transcript]
            prot += [f">{g.nuc.replace('NM', 'NP')} test protein", g.protein]
            if not g.starts:
                continue
            chrom = recipe._PRIMARY[g.nc.split(".")[0]]
            tag = f'tag "{g.status}";'
            attrs = f'gene_id "G{n}"; transcript_id "{g.nuc}"; {tag}'
            spans = self.spans(g)
            for start, end in spans:
                gtf.append(
                    f"chr{chrom}\tBestRefSeq\texon\t{start}\t{end}\t.\t{g.strand}\t.\t{attrs}"
                )
            lo = len(g.utr5) + 1  # first CDS base, in the transcript
            hi = len(g.utr5) + len(g.cds) - 3  # last CDS base: the stop codon is not CDS
            at = 1
            for start, end in spans:
                t0, t1 = at, at + end - start
                at = t1 + 1
                a, b = max(lo, t0), min(hi, t1)
                if a > b:
                    continue
                if g.strand == "+":
                    cs, ce = start + a - t0, start + b - t0
                else:
                    cs, ce = end - (b - t0), end - (a - t0)
                gtf.append(f"chr{chrom}\tBestRefSeq\tCDS\t{cs}\t{ce}\t.\t{g.strand}\t0\t{attrs}")
            gtf += g.extra_gtf
        genome = []
        for nc, bases in self.chrom.items():
            genome.append(f">{nc} Homo sapiens chromosome, GRCh38.p14 Primary Assembly")
            text = "".join(bases)
            genome += [text[i : i + 60] for i in range(0, len(text), 60)]
        genome.insert(0, ">NT_000000.1 an unrelated contig\nACGT")
        files = {
            recipe.SUMMARY: "\n".join(summary) + "\n",
            recipe.GTF: "\n".join(gtf) + "\n",
            recipe.RNA: "\n".join(rna) + "\n",
            recipe.PROTEIN: "\n".join(prot) + "\n",
            recipe.GENOME_NAME: "\n".join(genome) + "\n",
        }
        out = {}
        for name, text in files.items():
            out[name] = d / name
            out[name].write_bytes(gzip.compress(text.encode(), mtime=0))
        return out


def rows(world: World) -> recipe.Built:
    f = world.write()
    return recipe.read_rows(
        f[recipe.SUMMARY], f[recipe.GTF], f[recipe.RNA], f[recipe.PROTEIN], f[recipe.GENOME_NAME]
    )


def column(row: tuple[str, ...], name: str, columns: tuple[str, ...]) -> str:
    return row[columns.index(name)]


def test_a_minus_strand_transcript_joins_its_windows_back_to_its_sequence(tmp_path: Path) -> None:
    """Exon 1 of a minus-strand gene holds the HIGHEST genomic coordinates, reverse-complemented."""
    world = World(tmp_path)
    gene = world.add(
        Gene(
            "HGNC:9",
            "NM_000001.1",
            "MKTARSL",
            strand="-",
            cuts=(16, 31),
            starts=(300, 700, 1500),  # genome order: exon 3 of the transcript comes first
        )
    )
    built = rows(world)
    assert built.dropped == []
    (t_row,) = built.transcripts
    assert column(t_row, "strand", recipe.TRANSCRIPT_COLUMNS) == "-"
    assert column(t_row, "chrom", recipe.TRANSCRIPT_COLUMNS) == "12"
    exons = built.exons
    assert [int(e[2]) for e in exons] == [1, 2, 3]
    lengths = gene.exon_lengths()
    genomic = [(int(e[5]), int(e[6])) for e in exons]
    assert genomic == world.spans(gene)
    assert genomic[0][0] > genomic[1][0] > genomic[2][0]
    joined = ""
    for e in exons:
        window, w_start = e[9], int(e[7])
        core = window[int(e[5]) - w_start : int(e[6]) - w_start + 1]
        joined += revcomp(core)
        assert window == window.upper(), "the soft-masked genome is upper-cased"
        assert len(window) == int(e[8]) - w_start + 1
        assert int(e[5]) - w_start == recipe.FLANK and int(e[8]) - int(e[6]) == recipe.FLANK
    assert joined == gene.transcript == t_row[9]
    assert [(int(e[3]), int(e[4])) for e in exons] == [
        (1, lengths[0]),
        (lengths[0] + 1, lengths[0] + lengths[1]),
        (lengths[0] + lengths[1] + 1, sum(lengths)),
    ]


def test_cds_bounds_include_the_stop_codon(tmp_path: Path) -> None:
    """cds_end is the stop codon's last base, also where an exon boundary splits the stop."""
    world = World(tmp_path)
    protein = "MKTARSL"
    utr5 = "GGAATTCCGG"  # 10 bases
    start_at = len(utr5) + 1
    stop_last = len(utr5) + 3 * len(protein) + 3
    # the first exon ends one base into the stop codon (TAA)
    cut = stop_last - 2
    world.add(Gene("HGNC:1", "NM_000001.1", protein, cuts=(cut,), starts=(200, 900)))
    # a minus-strand twin with the cut in the same place
    world.add(Gene("HGNC:2", "NM_000002.1", protein, strand="-", cuts=(cut,), starts=(1500, 1800)))
    built = rows(world)
    assert built.dropped == []
    for t_row in built.transcripts:
        assert int(t_row[7]) == start_at
        assert int(t_row[8]) == stop_last == int(t_row[7]) + 3 * len(protein) + 2
        assert t_row[9][int(t_row[8]) - 3 : int(t_row[8])] == "TAA"


def test_selenocysteine_and_a_non_atg_start_pass_and_a_frameshift_is_dropped_with_its_reason(
    tmp_path: Path,
) -> None:
    world = World(tmp_path)
    # U against TGA and a non-ATG start for an M are allowed (the atlas's MAN002 does too)
    world.add(Gene("HGNC:1", "NM_000001.1", "MKUAR", starts=(100,)))
    world.add(Gene("HGNC:2", "NM_000002.1", "MKTAR", coding="CTG" + "AAAACCGCCCGG", starts=(400,)))
    # a base missing from the coding sequence shifts every later codon: the frame is lost
    world.add(Gene("HGNC:3", "NM_000003.1", "MKTAR", coding="ATGAAAACCGCCCG", starts=(700,)))
    # U against a codon that is not TGA is NOT allowed, nor is any other residue mismatch
    world.add(Gene("HGNC:4", "NM_000004.1", "MUTAR", coding="ATGAAAACCGCCCGG", starts=(1000,)))
    world.add(Gene("HGNC:5", "NM_000005.1", "MKTAR", coding="ATGCTGACCGCCCGG", starts=(1300,)))
    built = rows(world)
    assert [t[0] for t in built.transcripts] == ["HGNC:1", "HGNC:2"]
    assert {e[0] for e in built.exons} == {"HGNC:1", "HGNC:2"}
    by_reason = {d["reason"]: d for d in built.dropped}
    assert set(by_reason) == {"MAN002: the coding sequence does not translate to the protein"}
    entry = built.dropped[0]
    assert entry["count"] == 3
    assert entry["examples"] == [
        "HGNC:3 NM_000003.1: cds length 17 is not a multiple of 3",
        "HGNC:4 NM_000004.1: residue 2 is U in the protein but AAA encodes K",
        "HGNC:5 NM_000005.1: residue 2 is K in the protein but CTG encodes L",
    ]


def test_a_window_near_a_chromosome_end_is_clipped_and_states_its_own_bounds(
    tmp_path: Path,
) -> None:
    world = World(tmp_path, length=400)
    gene = world.add(Gene("HGNC:1", "NM_000001.1", "MKTAR", cuts=(20,), starts=(10, 360)))
    built = rows(world)
    assert built.dropped == []
    first, second = built.exons
    n = recipe.FLANK
    assert (first[5], first[7], first[8]) == ("10", "1", str(int(first[6]) + n))
    assert len(first[9]) == int(first[8]) - 1 + 1
    last_end = world.spans(gene)[1][1]
    assert (second[8], second[7]) == ("400", str(360 - n))
    assert last_end <= 400 and len(second[9]) == 400 - (360 - n) + 1


def test_an_ambiguity_code_in_a_flank_becomes_n_and_in_an_exon_core_drops_the_transcript(
    tmp_path: Path,
) -> None:
    world = World(tmp_path)
    ok = world.add(Gene("HGNC:1", "NM_000001.1", "MKTAR", starts=(300,)))
    bad = world.add(Gene("HGNC:2", "NM_000002.1", "MKTAR", starts=(900,)))
    (s1, e1), (s2, _) = world.spans(ok)[0], world.spans(bad)[0]
    r_at, n_at, k_at = s1 - 3, e1 + 4, e1 + 6  # 1-based positions, all in the flanks
    world.chrom["NC_000012.12"][r_at - 1] = "R"
    world.chrom["NC_000012.12"][k_at - 1] = "K"  # one in each flank: both are rewritten
    world.chrom["NC_000012.12"][n_at - 1] = "n"  # an N already there: counted separately
    world.chrom["NC_000012.12"][s2 + 1] = "Y"  # the third base of the second gene's core
    built = rows(world)
    assert [t[0] for t in built.transcripts] == ["HGNC:1"]
    (exon,) = built.exons
    window, w_start = exon[9], int(exon[7])
    assert {window[p - w_start] for p in (r_at, n_at, k_at)} == {"N"}
    assert set(window) <= set("ACGTN")
    assert (built.flank_iupac, built.flank_n) == (2, 1)
    assert built.dropped == [
        {
            "reason": "MAN001: the exon cores do not reproduce the RNA sequence",
            "count": 1,
            "examples": ["HGNC:2 NM_000002.1: an exon core holds a base other than A, C, G or T"],
        }
    ]
    assert bad.hgnc == "HGNC:2"


def test_only_the_mane_select_transcript_of_a_gene_is_read(tmp_path: Path) -> None:
    world = World(tmp_path)
    gene = world.add(Gene("HGNC:1", "NM_000001.1", "MKTAR", starts=(300,)))
    (start, end) = world.spans(gene)[0]
    # a Plus Clinical row in the summary, and a stray exon of the Select transcript's id that
    # the GTF tags Plus Clinical: neither may be used
    world.add(Gene("HGNC:1", "NM_000009.1", "MKTAR", status="MANE Plus Clinical"))
    stray = (
        f"chr12\tBestRefSeq\texon\t{end + 100}\t{end + 120}\t.\t+\t.\t"
        f'gene_id "G1"; transcript_id "NM_000001.1"; tag "MANE Plus Clinical";'
    )
    gene.extra_gtf.append(stray)
    built = rows(world)
    assert built.dropped == []
    assert [t[1] for t in built.transcripts] == ["NM_000001.1"]
    assert len(built.exons) == 1 and built.rows_read == 1
    assert (built.exons[0][5], built.exons[0][6]) == (str(start), str(end))


def test_a_non_primary_accession_is_listed_and_does_not_count_toward_the_ceiling(
    tmp_path: Path,
) -> None:
    world = World(tmp_path)
    world.add(Gene("HGNC:1", "NM_000001.1", "MKTAR", starts=(300,)))
    for n in range(2, 2 + recipe.MAX_FAILURES + 3):
        world.add(Gene(f"HGNC:{n}", f"NM_{n:06d}.1", "MKTAR", nc="NT_187633.1"))
    built = rows(world)
    assert [t[0] for t in built.transcripts] == ["HGNC:1"]
    (entry,) = built.dropped
    assert entry["count"] == recipe.MAX_FAILURES + 3
    assert entry["reason"].startswith("MANE Select transcript placed on a non-primary accession")
    assert entry["examples"][0] == "HGNC:10 NM_000010.1 on NT_187633.1"
    assert built.rows_read == recipe.MAX_FAILURES + 4


def test_a_non_coding_transcript_and_a_gene_without_an_hgnc_id_are_listed_not_counted(
    tmp_path: Path,
) -> None:
    world = World(tmp_path)
    world.add(Gene("HGNC:1", "NM_000001.1", "MKTAR", starts=(300,)))
    for n in range(recipe.MAX_FAILURES + 2):
        world.add(Gene(f"HGNC:{100 + n}", f"NR_{n:06d}.1", "MKTAR"))
        world.add(Gene("", f"NM_{900000 + n}.1", "MKTAR"))
    built = rows(world)
    assert [t[0] for t in built.transcripts] == ["HGNC:1"]
    counts = {d["reason"].split(",")[0].split(" (")[0]: d["count"] for d in built.dropped}
    assert counts == {
        "MANE Select non-coding transcript": 52,
        "MANE Select transcript of a gene with no HGNC id": 52,
    }
    assert {d["examples"][0] for d in built.dropped} == {
        "HGNC:100 NR_000000.1",
        "NM_900000.1 (blank)",
    }
    assert recipe._PRIMARY["NC_000023"] == "X" and recipe._PRIMARY["NC_012920"] == "MT"


def _failing(world: World, count: int) -> None:
    for n in range(2, 2 + count):
        world.add(
            Gene(f"HGNC:{n}", f"NM_{n:06d}.1", "MKTAR", coding="ATGAAAACCGCCCG", starts=(100 * n,))
        )


def test_more_than_fifty_failures_refuse_and_fifty_do_not(tmp_path: Path) -> None:
    ok = World(tmp_path / "a", length=9000)
    (tmp_path / "a").mkdir()
    ok.tmp = tmp_path / "a"
    ok.add(Gene("HGNC:1", "NM_000001.1", "MKTAR", starts=(50,)))
    _failing(ok, recipe.MAX_FAILURES)
    built = rows(ok)
    assert built.dropped[0]["count"] == recipe.MAX_FAILURES
    (tmp_path / "b").mkdir()
    over = World(tmp_path / "b", length=9000)
    over.add(Gene("HGNC:1", "NM_000001.1", "MKTAR", starts=(50,)))
    _failing(over, recipe.MAX_FAILURES + 1)
    with pytest.raises(PackageRefused, match="more than 50"):
        rows(over)


@pytest.fixture
def pinned(monkeypatch: pytest.MonkeyPatch) -> object:
    def _pin(files: dict[str, Path]) -> None:
        monkeypatch.setattr(
            recipe,
            "PINNED_SHA256",
            {n: hashlib.sha256(p.read_bytes()).hexdigest() for n, p in files.items()},
        )

    return _pin


def _build(files: dict[str, Path], out: Path) -> Path:
    return recipe.build(
        files[recipe.SUMMARY],
        files[recipe.GTF],
        files[recipe.RNA],
        files[recipe.PROTEIN],
        files[recipe.GENOME_NAME],
        retrieved_on=date(2026, 11, 17),
        revision=1,
        built_by=BUILT_BY,
        out=out,
        commit_of=lambda: RECIPE_COMMIT,
    )


def _two_genes(tmp_path: Path) -> World:
    world = World(tmp_path)
    world.add(Gene("HGNC:9", "NM_000002.1", "MKTAR", strand="-", cuts=(20,), starts=(300, 900)))
    world.add(Gene("HGNC:10", "NM_000001.1", "MKTARSL", nc="NC_000023.11", starts=(500,)))
    return world


def test_the_package_states_its_inputs_and_orders_rows_by_gene_then_exon(
    tmp_path: Path, pinned: object
) -> None:
    files = _two_genes(tmp_path).write()
    pinned(files)  # type: ignore[operator]
    out = _build(files, tmp_path / "pkg")
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert (manifest["package"], manifest["version"]) == ("mane_select", "1.5")
    assert (manifest["ids"], manifest["coverage"]) == ("hgnc_id", "every_gene")
    assert manifest["input_counts"] == {"rows": 2} and manifest["dropped"] == []
    assert [(t["name"], t["rows"], t["genes"]) for t in manifest["tables"]] == [
        ("mane_exons", 3, 2),
        ("mane_transcripts", 2, 2),
    ]
    names = [i["name"] for i in manifest["inputs"]]
    assert names == [recipe.SUMMARY, recipe.GTF, recipe.RNA, recipe.PROTEIN, recipe.GENOME_NAME]
    assert {i["retrieved_on"] for i in manifest["inputs"]} == {"2026-11-17"}
    assert manifest["inputs"][0]["url"].endswith("release_1.5/" + recipe.SUMMARY)
    assert manifest["source"]["version"] == "MANE v1.5"
    raw = gzip.decompress((out / "mane_transcripts.tsv.gz").read_bytes()).decode()
    lines = [line.split("\t") for line in raw.splitlines()]
    assert lines[0] == list(recipe.TRANSCRIPT_COLUMNS)
    assert [(line[0], line[5], line[11]) for line in lines[1:]] == [
        ("HGNC:10", "X", "mane_select"),  # NC_000023 is X; "HGNC:10" sorts before "HGNC:9"
        ("HGNC:9", "12", "mane_select"),
    ]
    raw_e = gzip.decompress((out / "mane_exons.tsv.gz").read_bytes()).decode()
    assert [(r.split("\t")[0], r.split("\t")[2]) for r in raw_e.splitlines()[1:]] == [
        ("HGNC:10", "1"),
        ("HGNC:9", "1"),
        ("HGNC:9", "2"),
    ]


def test_other_bytes_are_refused_and_nothing_is_written(tmp_path: Path) -> None:
    files = _two_genes(tmp_path).write()
    with pytest.raises(PackageRefused, match="not the pinned"):
        _build(files, tmp_path / "pkg")
    assert not (tmp_path / "pkg").exists()


def test_one_changed_input_is_refused_by_name(tmp_path: Path, pinned: object) -> None:
    files = _two_genes(tmp_path).write()
    pinned(files)  # type: ignore[operator]
    files[recipe.PROTEIN].write_bytes(gzip.compress(b">NP_000001.1\nMK\n", mtime=0))
    with pytest.raises(PackageRefused, match=recipe.PROTEIN):
        _build(files, tmp_path / "pkg")
    assert not (tmp_path / "pkg").exists()


def test_the_pins_cover_exactly_the_five_inputs() -> None:
    assert (
        set(recipe.PINNED_SHA256)
        == set(recipe.DOWNLOAD_URLS)
        == {
            recipe.SUMMARY,
            recipe.GTF,
            recipe.RNA,
            recipe.PROTEIN,
            recipe.GENOME_NAME,
        }
    )
    assert all(len(v) == 64 for v in recipe.PINNED_SHA256.values())
    assert recipe.FLANK == 50 and recipe.UPSTREAM_RELEASE == "1.5"


def test_the_real_v1_5_files_build_when_present(tmp_path: Path) -> None:
    """The real build streams a 3 GB genome (minutes), so it also needs FEEDS_REAL=1."""
    local = Path(__file__).resolve().parents[3] / "local" / "mane"
    names = (recipe.SUMMARY, recipe.GTF, recipe.RNA, recipe.PROTEIN, recipe.GENOME_NAME)
    if not all((local / n).is_file() for n in names):
        pytest.skip(f"{local} does not hold the MANE v1.5 files (gitignored; a curator's copy)")
    if os.environ.get("FEEDS_REAL") != "1":
        pytest.skip("set FEEDS_REAL=1 to build from the real files (takes minutes)")
    built = recipe.read_rows(*(local / n for n in names))
    dropped = sum(d["count"] for d in built.dropped)
    assert built.rows_read == 19363
    assert len(built.transcripts) == 19363 - dropped == 19179
    assert {d["count"] for d in built.dropped} == {1, 61, 52, 70}
    by_gene = {t[0]: t for t in built.transcripts}
    cols = recipe.TRANSCRIPT_COLUMNS
    tbx5 = next(
        t for t in built.transcripts if column(t, "refseq_transcript", cols) == "NM_181486.4"
    )
    assert column(tbx5, "strand", cols) == "-" and column(tbx5, "cds_start", cols) == "556"
    assert column(tbx5, "protein_sequence", cols)[278] == "R"
    nr2f2 = next(
        t for t in by_gene.values() if column(t, "refseq_transcript", cols) == "NM_021005.4"
    )
    assert column(nr2f2, "cds_start", cols) == "1407"
    assert json.dumps(built.dropped)  # serialisable
