"""Weekly drift check: has an upstream source changed since its baseline was captured?

Usage: drift.py [--baselines DIR] [--open-issues | --capture]

For each source in ``SOURCES`` this runs a probe, compares the fingerprint it returns
with ``sources/<id>/drift.json`` and prints a report. Exit codes: 0 every source
matches its baseline, 1 at least one drifted, 2 at least one probe failed or has no
baseline. A failure outranks a drift, and is never reported as clean: a probe that
could not run says nothing about its source.

A baseline means "nothing changed since <captured_on>". It does not say the pinned
package is current; the package may predate the baseline (ClinGen and GenCC do).

Four sources use hvantk's ``drift`` probe, read-only, through its ``--json`` output
(installed at a pinned commit by the workflow); the two Orphanet XML products use
``orphadata_probe``, which reads the ``<JDBOR>`` version and date from the first 4 KB.
hvantk's own verdict, ``expected`` and ``diff`` are ignored: its ``expected`` is a
fingerprint packaged inside hvantk, not this repository's baseline.
"""

from __future__ import annotations

import argparse
import datetime
import json
import re
import subprocess
import sys
import urllib.request
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
REPO = "chd-atlas/feeds"
HEADER_BYTES = 4096
HVANTK_TIMEOUT = 180  # seconds, passed to ``hvantk drift --timeout``
# Keys of hvantk's ``observed`` that identify the upstream bytes. ``fetched_at`` and
# ``informational`` (e.g. Last-Modified) change on every run and are never compared.
COMPARED = ("probe_version", "source_version", "headers", "checksums", "extras")

Fetch = Callable[[str, int], bytes]
Hvantk = Callable[[str, int], tuple[int, str]]
Gh = Callable[[list[str]], str]


class ProbeFailed(Exception):
    """A probe could not produce a fingerprint."""


@dataclass(frozen=True)
class HvantkSource:
    dataset: str
    package: str  # the pinned package this baseline is read against, ``<id>/<version>``


@dataclass(frozen=True)
class OrphadataSource:
    url: str
    package: str


Source = HvantkSource | OrphadataSource

SOURCES: dict[str, Source] = {
    "hgnc": HvantkSource("hgnc:lookup", "hgnc/2026-10-06"),
    "clingen": HvantkSource("clingen:gene-disease", "clingen/2026-08-03"),
    "gencc": HvantkSource("gencc:submissions", "gencc/88c77272"),
    "expression_atlas_e_mtab_6814": HvantkSource(
        "expression-atlas:dataset", "expression_atlas_e_mtab_6814/92484a97"
    ),
    "orphanet_product4": OrphadataSource(
        "https://www.orphadata.com/data/xml/en_product4.xml", "orphanet_product4/4f44e8a6"
    ),
    "orphanet_product1": OrphadataSource(
        "https://www.orphadata.com/data/xml/en_product1.xml", "orphanet_product1/df8d562a"
    ),
}


def compare(baseline: Mapping[str, Any], observed: Mapping[str, Any]) -> list[str]:
    """The keys whose values differ between two fingerprints, sorted.

    Nested mappings are walked and reported as ``outer/inner`` (``checksums/f.csv``).
    A key on one side only is a change. ``fetched_at`` and ``informational`` are skipped
    at the top level whichever side carries them; for hvantk fingerprints every other
    key (``COMPARED``) is compared.
    """
    left, right = _flat(baseline), _flat(observed)
    return sorted(
        k for k in set(left) | set(right) if k not in right or k not in left or left[k] != right[k]
    )


def _flat(value: Mapping[str, Any], prefix: str = "") -> dict[str, Any]:
    """Leaf values keyed ``outer/inner``; top-level ``fetched_at`` and ``informational`` dropped."""
    flat: dict[str, Any] = {}
    for key, item in value.items():
        if not prefix and key in ("fetched_at", "informational"):
            continue
        if isinstance(item, Mapping) and item:
            flat |= _flat(item, f"{prefix}{key}/")
        else:
            flat[f"{prefix}{key}"] = item
    return flat


def fingerprint_of_hvantk(observed: Mapping[str, Any]) -> dict[str, Any]:
    """The part of hvantk's ``observed`` that names the upstream bytes."""
    return {k: observed[k] for k in COMPARED if k in observed}


def orphadata_probe(url: str, fetch: Fetch) -> dict[str, str]:
    """The ``version`` and ``date`` on the root ``<JDBOR>`` element, from the first 4 KB.

    The prefix is cut mid-document, so it is matched with a regular expression rather than
    parsed. A reply that is not a ``<JDBOR>`` document, or states either value blank, is a
    failed probe: recording a blank release would let a changed file compare equal.
    """
    try:
        head = fetch(url, HEADER_BYTES)
    except Exception as exc:
        raise ProbeFailed(f"{url}: {type(exc).__name__}: {exc}") from exc
    text = head.decode("utf-8", errors="replace")
    root = re.search(r"<JDBOR\b([^>]*)>", text)
    if root is None:
        raise ProbeFailed(f"{url}: no <JDBOR> element in the first {HEADER_BYTES} bytes")
    attrs = dict(re.findall(r'([\w:-]+)="([^"]*)"', root.group(1)))
    version, date = attrs.get("version", "").strip(), attrs.get("date", "").strip()
    if not version or not date:
        raise ProbeFailed(f"{url}: <JDBOR> states no version or no date")
    return {"version": version, "date": date}


def http_fetch(url: str, nbytes: int) -> bytes:
    """The first ``nbytes`` of a URL (a Range request; reads at most ``nbytes`` regardless)."""
    request = urllib.request.Request(url, headers={"Range": f"bytes=0-{nbytes - 1}"})
    with urllib.request.urlopen(request, timeout=60) as response:  # noqa: S310
        return bytes(response.read(nbytes))


def run_hvantk(dataset: str, timeout: int) -> tuple[int, str]:
    done = subprocess.run(
        ["hvantk", "drift", dataset, "--json", "--timeout", str(timeout)],
        capture_output=True,
        text=True,
        timeout=timeout + 120,
        check=False,
    )
    return done.returncode, done.stdout


def run_gh(argv: list[str]) -> str:
    done = subprocess.run(["gh", *argv], capture_output=True, text=True, check=True)
    return done.stdout


def _hvantk_probe(dataset: str, hvantk: Hvantk) -> dict[str, Any]:
    try:
        code, out = hvantk(dataset, HVANTK_TIMEOUT)
    except Exception as exc:
        raise ProbeFailed(f"hvantk drift {dataset}: {type(exc).__name__}: {exc}") from exc
    try:
        entries = json.loads(out)
        entry = entries[0]
    except (ValueError, IndexError, KeyError, TypeError) as exc:
        raise ProbeFailed(f"hvantk drift {dataset}: exit {code}, unreadable output") from exc
    if not isinstance(entry, dict):
        raise ProbeFailed(f"hvantk drift {dataset}: exit {code}, unreadable output")
    if entry.get("status") == "probe_failed":
        raise ProbeFailed(f"hvantk drift {dataset}: {entry.get('probe_error') or 'probe_failed'}")
    observed = entry.get("observed")
    if code not in (0, 1) or not isinstance(observed, dict) or not observed:
        raise ProbeFailed(f"hvantk drift {dataset}: exit {code} with no observed fingerprint")
    return fingerprint_of_hvantk(observed)


def probe(spec: Source, hvantk: Hvantk, fetch: Fetch) -> dict[str, Any]:
    if isinstance(spec, HvantkSource):
        return _hvantk_probe(spec.dataset, hvantk)
    return dict(orphadata_probe(spec.url, fetch))


def _issue_body(source_id: str, baseline: Mapping[str, Any], keys: list[str]) -> str:
    listed = "\n".join(f"- `{k}`" for k in keys)
    return (
        f"The upstream source of `{source_id}` changed since its baseline was captured.\n\n"
        f"- Pinned package: `{baseline.get('package', 'unknown')}`\n"
        f"- Baseline captured on: {baseline.get('captured_on', 'unknown')} "
        "(nothing changed since that date was the baseline's claim; it does not say the "
        "pinned package is current)\n"
        f"- Changed fingerprint keys:\n{listed}\n\n"
        "Decide whether to rebuild the package, then refresh the baseline with "
        "`scripts/drift.py --capture`.\n"
    )


def _open_issue(gh: Gh, source_id: str, body: str) -> None:
    title = f"drift: {source_id}"
    found = json.loads(
        gh(
            ["issue", "list", "--repo", REPO, "--state", "open", "--search", f"{title} in:title"]
            + ["--json", "number,title"]
        )
    )
    numbers = [str(i["number"]) for i in found if i.get("title") == title]
    if numbers:
        gh(["issue", "edit", numbers[0], "--repo", REPO, "--body", body])
    else:
        gh(["issue", "create", "--repo", REPO, "--title", title, "--body", body])


def main(
    argv: Sequence[str] | None = None,
    *,
    hvantk: Hvantk = run_hvantk,
    fetch: Fetch = http_fetch,
    gh: Gh = run_gh,
    sources: Mapping[str, Source] | None = None,
    today: str | None = None,
) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0] if __doc__ else None)
    parser.add_argument("--baselines", type=Path, default=ROOT / "sources")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--open-issues", action="store_true")
    group.add_argument("--capture", action="store_true", help="write every baseline, today")
    args = parser.parse_args(argv)
    specs = SOURCES if sources is None else sources

    observed: dict[str, dict[str, Any]] = {}
    failures: dict[str, str] = {}
    for source_id in sorted(specs):
        try:
            observed[source_id] = probe(specs[source_id], hvantk, fetch)
        except ProbeFailed as exc:
            failures[source_id] = str(exc)

    if args.capture:
        for source_id, message in sorted(failures.items()):
            print(f"{source_id}: PROBE FAILED: {message}")
        if failures:
            print("nothing written: a baseline is captured only when every probe answered")
            return 2
        stamp = today or datetime.date.today().isoformat()
        for source_id, fp in sorted(observed.items()):
            record = {"captured_on": stamp, "package": specs[source_id].package, "fingerprint": fp}
            path = args.baselines / source_id / "drift.json"
            path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            print(f"{source_id}: baseline captured {stamp}")
        return 0

    drifted = 0
    for source_id in sorted(specs):
        if source_id in failures:
            print(f"{source_id}: PROBE FAILED: {failures[source_id]}")
            continue
        path = args.baselines / source_id / "drift.json"
        try:
            baseline = json.loads(path.read_text())
            keys = compare(baseline["fingerprint"], observed[source_id])
        except (OSError, ValueError, KeyError, TypeError) as exc:
            failures[source_id] = f"no baseline: {path} unreadable ({type(exc).__name__})"
            print(f"{source_id}: PROBE FAILED: {failures[source_id]}")
            continue
        if not keys:
            print(f"{source_id}: clean (nothing changed since {baseline.get('captured_on')})")
            continue
        drifted += 1
        print(f"{source_id}: DRIFTED since {baseline.get('captured_on')}: {', '.join(keys)}")
        if args.open_issues:
            _open_issue(gh, source_id, _issue_body(source_id, baseline, keys))
    if failures:
        return 2
    return 1 if drifted else 0


if __name__ == "__main__":
    sys.exit(main())
