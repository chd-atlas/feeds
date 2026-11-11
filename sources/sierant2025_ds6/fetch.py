"""Download the Sierant 2025 Dataset S6 workbook from Europe PMC, standard library only.

pnas.org answers 403 to a scripted request (as recorded in the atlas's
``mirrors/sources.yaml``), so the workbook is taken from the Europe PMC
open-access supplementary bundle for PMC12002227, and its sha256 is what proves
it is the file the atlas pinned.

A downloaded bundle is untrusted: only the one named member is read, its
bytes are written under a name this module chooses, and the sha256 is checked
before the file is left in ``dest``.
"""

from __future__ import annotations

import hashlib
import shutil
import tempfile
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath

BUNDLE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12002227/supplementaryFiles"
WORKBOOK_NAME = "pnas.2420343122.sd06.xlsx"
PINNED_SHA256 = "d8e9f7777b329c810adde8b9f4a0f91a4d1a6211b760daeb9329ae6f36606daa"
# A declared uncompressed size above this refuses the member before it is read. The real
# workbook is 233,828 bytes (measured 2026-11-11), so this leaves a factor of about 21.
# zipfile reads no more than the declared size, so a header that understates the size
# cannot get past the cap and then inflate beyond it.
MAX_MEMBER_BYTES = 5_000_000
# Seconds without a byte before the download is abandoned.
TIMEOUT = 120


class FetchRefused(Exception):
    """The bundle is not the one pinned (or is unsafe); nothing is left in ``dest``."""


def fetch(dest: Path) -> Path:
    """Download the bundle, extract only the workbook into ``dest`` and return its path."""
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / WORKBOOK_NAME
    if target.exists():
        raise FetchRefused(f"{target} exists; remove it first")
    with tempfile.TemporaryDirectory(prefix=".bundle.") as scratch:
        bundle = Path(scratch) / "bundle.zip"
        with (
            urllib.request.urlopen(BUNDLE_URL, timeout=TIMEOUT) as response,
            bundle.open("wb") as handle,
        ):
            shutil.copyfileobj(response, handle)
        try:
            archive = zipfile.ZipFile(bundle)
        except zipfile.BadZipFile as exc:
            raise FetchRefused(f"{BUNDLE_URL} did not return a zip archive: {exc}") from exc
        with archive:
            matches = []
            for info in archive.infolist():
                path = PurePosixPath(info.filename)
                if path.name != WORKBOOK_NAME:
                    continue
                if path.is_absolute() or ".." in path.parts or "\\" in info.filename:
                    raise FetchRefused(
                        f"bundle member '{info.filename}' has a path that escapes the "
                        "destination; refusing the bundle"
                    )
                if info.file_size > MAX_MEMBER_BYTES:
                    raise FetchRefused(
                        f"bundle member '{info.filename}' declares {info.file_size} bytes, "
                        f"over the {MAX_MEMBER_BYTES}-byte cap"
                    )
                matches.append(info)
            if len(matches) != 1:
                raise FetchRefused(
                    f"expected exactly one member named {WORKBOOK_NAME} in the bundle, "
                    f"found {len(matches)}"
                )
            data = archive.read(matches[0])
    digest = hashlib.sha256(data).hexdigest()
    if digest != PINNED_SHA256:
        raise FetchRefused(
            f"{WORKBOOK_NAME} from the bundle has sha256 {digest}, not the pinned {PINNED_SHA256}; "
            "nothing was written"
        )
    target.write_bytes(data)
    return target
