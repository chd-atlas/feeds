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
from collections.abc import Callable
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
    in_rna: bool = True
    in_protein: bool = True
    ens_nuc: str | None = None  # the summary's Ensembl transcript, when not well-formed by default
    gtf_chrom: str | None = None  # the GTF seqname (without "chr"), when not the summary's
    gtf_strand: str | None = None

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
        self.omit_from_genome: set[str] = set()
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
                f"\t{g.ens_nuc or f'ENST{num}.1'}\tENSP{num}.1\t{g.status}\t{g.nc}\t1\t2"
                f"\t{g.strand}"
            )
            if g.in_rna:
                rna += [f">{g.nuc} Homo sapiens test, mRNA", g.transcript]
            if g.in_protein:
                prot += [f">{g.nuc.replace('NM', 'NP')} test protein", g.protein]
            if not g.starts:
                continue
            chrom = g.gtf_chrom or recipe._PRIMARY[g.nc.split(".")[0]]
            strand = g.gtf_strand or g.strand
            tag = f'tag "{g.status}";'
            attrs = f'gene_id "G{n}"; transcript_id "{g.nuc}"; {tag}'
            spans = self.spans(g)
            for start, end in spans:
                gtf.append(f"chr{chrom}\tBestRefSeq\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs}")
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
                gtf.append(f"chr{chrom}\tBestRefSeq\tCDS\t{cs}\t{ce}\t.\t{strand}\t0\t{attrs}")
            gtf += g.extra_gtf
        genome = []
        for nc, bases in self.chrom.items():
            if nc in self.omit_from_genome:
                continue
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
    assert entry["reason"].startswith("MANE Select transcript placed only on an alternate-locus")
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


def _good(world: World) -> None:
    world.add(Gene("HGNC:1", "NM_000001.1", "MKTAR", starts=(300,)))


def _bad(n: int, **kw: object) -> Gene:
    """The n-th misbehaving gene: its own id, accession and place on the chromosome."""
    fields: dict[str, object] = {"starts": (1000 + 150 * n,)}
    fields.update(kw)
    return Gene(f"HGNC:{n + 2}", f"NM_{n + 2:06d}.1", "MKTAR", **fields)  # type: ignore[arg-type]


def _case_chrom_absent(world: World, n: int) -> None:
    world.add(_bad(n, nc="NC_000013.11"))
    world.omit_from_genome.add("NC_000013.11")


def _case_no_rna(world: World, n: int) -> None:
    world.add(_bad(n, in_rna=False))


def _case_no_protein(world: World, n: int) -> None:
    world.add(_bad(n, in_protein=False))


def _case_bad_ensembl(world: World, n: int) -> None:
    world.add(_bad(n, ens_nuc="ENST1.1"))


def _case_gtf_wrong_chromosome(world: World, n: int) -> None:
    world.add(_bad(n, gtf_chrom="13"))


def _case_gtf_wrong_strand(world: World, n: int) -> None:
    world.add(_bad(n, gtf_strand="-"))


def _case_no_exons(world: World, n: int) -> None:
    world.add(_bad(n, starts=()))


def _case_core_n(world: World, n: int) -> None:
    gene = world.add(_bad(n))
    world.chrom["NC_000012.12"][world.spans(gene)[0][0] + 1] = "N"


def _case_odd_accession(world: World, n: int) -> None:
    gene = _bad(n, starts=())
    gene.nuc = f"XM_{n + 2:06d}.1"
    world.add(gene)


def _case_odd_hgnc(world: World, n: int) -> None:
    gene = _bad(n, starts=())
    gene.hgnc = "HGNC:x"
    world.add(gene)


def _case_odd_chromosome(world: World, n: int) -> None:
    world.add(_bad(n, starts=(), nc="12"))


_NO_EXON = "no exon in the GTF on the summary's chromosome and strand"
# (case, reason, example for the first bad gene: HGNC:2, NM_000002.1)
_DROPS = [
    (
        _case_chrom_absent,
        "chromosome NC_000013.11 is not in the genome FASTA",
        "HGNC:2 NM_000002.1",
    ),
    (_case_no_rna, "no RNA sequence in the MANE FASTA", "HGNC:2 NM_000002.1"),
    (_case_no_protein, "no protein sequence in the MANE FASTA", "HGNC:2 NM_000002.1"),
    (
        _case_bad_ensembl,
        "ensembl_transcript does not fit the table's pattern",
        "HGNC:2 NM_000002.1",
    ),
    (_case_gtf_wrong_chromosome, _NO_EXON, "HGNC:2 NM_000002.1"),
    (_case_gtf_wrong_strand, _NO_EXON, "HGNC:2 NM_000002.1"),
    (_case_no_exons, _NO_EXON, "HGNC:2 NM_000002.1"),
    (
        _case_core_n,
        "MAN001: the exon cores do not reproduce the RNA sequence",
        "HGNC:2 NM_000002.1: an exon core holds a base other than A, C, G or T",
    ),
    (
        _case_odd_accession,
        "MANE Select transcript that is neither NM_ nor NR_",
        "HGNC:2 XM_000002.1",
    ),
    (
        _case_odd_hgnc,
        "MANE Select transcript whose HGNC id is not HGNC:<digits>",
        "NM_000002.1 (HGNC:x)",
    ),
    (
        _case_odd_chromosome,
        "MANE Select transcript on an unrecognised chromosome accession",
        "HGNC:2 NM_000002.1 on 12",
    ),
]
_IDS = [c[0].__name__.removeprefix("_case_") for c in _DROPS]


@pytest.mark.parametrize(("case", "reason", "example"), _DROPS, ids=_IDS)
def test_a_transcript_that_cannot_be_framed_is_dropped_alone_with_its_exact_reason(
    tmp_path: Path, case: object, reason: str, example: str
) -> None:
    """One bad transcript beside a good one: it is listed (never refused) and the other stays."""
    world = World(tmp_path, length=20000)
    _good(world)
    case(world, 0)  # type: ignore[operator]
    built = rows(world)
    assert [t[0] for t in built.transcripts] == ["HGNC:1"]
    assert built.dropped == [{"reason": reason, "count": 1, "examples": [example]}]


@pytest.mark.parametrize("case", [c[0] for c in _DROPS], ids=_IDS)
def test_every_such_drop_counts_toward_the_fifty_failure_ceiling(
    tmp_path: Path, case: object
) -> None:
    """51 of one kind refuse; 50 do not. Only the by-design exclusions are exempt."""
    for count in (recipe.MAX_FAILURES, recipe.MAX_FAILURES + 1):
        sub = tmp_path / str(count)
        sub.mkdir()
        world = World(sub, length=20000)
        _good(world)
        for n in range(count):
            case(world, n)  # type: ignore[operator]
        if count > recipe.MAX_FAILURES:
            with pytest.raises(PackageRefused, match="more than 50"):
                rows(world)
        else:
            assert rows(world).dropped[0]["count"] == count


def _edit_summary(world: World, edit: Callable[[list[str]], list[str]]) -> dict[str, Path]:
    files = world.write()
    path = files[recipe.SUMMARY]
    lines = gzip.decompress(path.read_bytes()).decode().splitlines()
    path.write_bytes(gzip.compress(("\n".join(edit(lines)) + "\n").encode(), mtime=0))
    return files


def _read_files(files: dict[str, Path]) -> recipe.Built:
    return recipe.read_rows(
        files[recipe.SUMMARY],
        files[recipe.GTF],
        files[recipe.RNA],
        files[recipe.PROTEIN],
        files[recipe.GENOME_NAME],
    )


def test_a_duplicate_hgnc_row_is_refused(tmp_path: Path) -> None:
    world = World(tmp_path)
    _good(world)
    files = _edit_summary(world, lambda lines: [*lines, lines[-1]])
    with pytest.raises(PackageRefused, match="HGNC:1 has more than one MANE Select row"):
        _read_files(files)


def test_a_summary_without_a_needed_column_is_refused(tmp_path: Path) -> None:
    world = World(tmp_path)
    _good(world)
    files = _edit_summary(
        world, lambda lines: [lines[0].replace("chr_strand", "strand"), *lines[1:]]
    )
    with pytest.raises(PackageRefused, match="no column chr_strand"):
        _read_files(files)


def test_a_short_summary_row_is_a_counted_failure_not_a_silent_skip(tmp_path: Path) -> None:
    world = World(tmp_path)
    _good(world)
    files = _edit_summary(world, lambda lines: [*lines, "GeneID:9\tENSG1\tHGNC:9"])
    built = _read_files(files)
    assert built.dropped == [
        {
            "reason": "summary row with too few cells to read its MANE status",
            "count": 1,
            "examples": ["line 3 (3 cells)"],
        }
    ]
    many = _edit_summary(world, lambda lines: [*lines, *["GeneID:9\tENSG1\tHGNC:9"] * 51])
    with pytest.raises(PackageRefused, match="more than 50"):
        _read_files(many)


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


def test_a_later_revision_reads_r2_and_the_genome_names_its_assembly(
    tmp_path: Path, pinned: object, capsys: pytest.CaptureFixture[str]
) -> None:
    world = _two_genes(tmp_path)
    world.chrom["NC_000012.12"][299 - 5] = "R"  # one flank IUPAC code for the operator's line
    files = world.write()
    pinned(files)  # type: ignore[operator]
    out = recipe.build(
        files[recipe.SUMMARY],
        files[recipe.GTF],
        files[recipe.RNA],
        files[recipe.PROTEIN],
        files[recipe.GENOME_NAME],
        retrieved_on=date(2026, 11, 17),
        revision=2,
        built_by=BUILT_BY,
        out=tmp_path / "pkg",
        commit_of=lambda: RECIPE_COMMIT,
    )
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["version"] == "1.5-r2"
    releases = {i["name"]: i["upstream_release"] for i in manifest["inputs"]}
    assert releases == {
        recipe.SUMMARY: "1.5",
        recipe.GTF: "1.5",
        recipe.RNA: "1.5",
        recipe.PROTEIN: "1.5",
        recipe.GENOME_NAME: "GCF_000001405.40 (GRCh38.p14)",
    }
    printed = capsys.readouterr().out
    assert "2 transcripts, 3 exons; 0 flank bases already N, 1 flank IUPAC codes" in printed


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
    assert {d["reason"]: d["count"] for d in built.dropped} == {
        "MAN002: the coding sequence does not translate to the protein": 1,
        "MANE Select non-coding transcript (NR_), which has no protein to frame": 70,
        "MANE Select transcript of a gene with no HGNC id (the table is keyed by HGNC id)": 52,
        recipe._NON_PRIMARY_REASON: 61,
    }
    assert len(built.exons) == 201933
    assert (built.flank_n, built.flank_iupac) == (1, 0)
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
    _codon_features_give_the_same_bounds(local / recipe.GTF, built)


def _codon_features_give_the_same_bounds(gtf: Path, built: recipe.Built) -> None:
    """Re-derive every kept transcript's bounds from start/stop_codon, joining split ones."""
    exons: dict[str, list[tuple[int, int, int, int]]] = {}
    for e in built.exons:
        exons.setdefault(e[1], []).append((int(e[5]), int(e[6]), int(e[3]), int(e[4])))
    strand = {t[1]: t[6] for t in built.transcripts}
    bounds = {t[1]: (int(t[7]), int(t[8])) for t in built.transcripts}
    found: dict[str, dict[str, list[int]]] = {}
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] not in ("start_codon", "stop_codon"):
                continue
            if 'tag "MANE Select";' not in cols[8]:
                continue
            name = cols[8].split('transcript_id "', 1)[1].split('"', 1)[0]
            if name in bounds:
                found.setdefault(name, {}).setdefault(cols[2], []).extend(
                    [int(cols[3]), int(cols[4])]
                )

    def at(name: str, pos: int) -> int:
        for g_start, g_end, t_start, _ in exons[name]:
            if g_start <= pos <= g_end:
                return t_start + (g_end - pos if strand[name] == "-" else pos - g_start)
        raise AssertionError(f"{name}: {pos} is in no exon")

    for name, (cds_start, cds_end) in bounds.items():
        start, stop = found[name]["start_codon"], found[name]["stop_codon"]
        first, last = (min(start), max(stop)) if strand[name] == "+" else (max(start), min(stop))
        assert (at(name, first), at(name, last)) == (cds_start, cds_end), name
