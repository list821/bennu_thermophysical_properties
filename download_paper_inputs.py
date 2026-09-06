"""Download public inputs used by the Bennu 2019 Nature Astronomy paper.

Default mode downloads the manageable core set: paper and instrument documents,
two OTES calibrated observations, their geometry, two representative OVIRS
products, a public 12.6 m Bennu shape model, and date-relevant SPICE kernels.

The complete OVIRS fit set contains 34,254 FITS products plus PDS labels and is
roughly 19 GB.  It is downloaded only when both --full-ovirs and
--confirm-large-download are supplied.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
USER_AGENT = "Bennu-ATPM-reproduction/1.0 (scientific data download)"

PAPER_URL = (
    "https://ntrs.nasa.gov/api/citations/20190025172/downloads/20190025172.pdf"
)
OVIRS = "https://sbnarchive.psi.edu/pds4/orex/orex.ovirs"
OTES = "https://sbnarchive.psi.edu/pds4/orex/orex.otes"
SPICE = "https://naif.jpl.nasa.gov/pub/naif/pds/pds4/orex/orex_spice/spice_kernels"


CORE_FILES: list[tuple[str, str, str]] = [
    (PAPER_URL, "paper/2019_bennu_properties_ntrs.pdf", "paper"),
    (f"{OVIRS}/document/ovirs_sis.pdf", "docs/ovirs_sis.pdf", "documentation"),
    (
        f"{OVIRS}/document/ovirs_cal_doc_v3.0.pdf",
        "docs/ovirs_cal_doc_v3.0.pdf",
        "documentation",
    ),
    (f"{OTES}/document/otes_sis.pdf", "docs/otes_sis.pdf", "documentation"),
    (
        f"{OTES}/document/otes_ops_timeline.pdf",
        "docs/otes_ops_timeline.pdf",
        "documentation",
    ),
]

for stem in (
    "20181108T040144S029_ote_scil2",
    "20181109T040143S209_ote_scil2",
):
    for suffix in ("dat", "xml"):
        CORE_FILES.append(
            (
                f"{OTES}/data_calibrated/approach/{stem}.{suffix}",
                f"otes/science/{stem}.{suffix}",
                "OTES calibrated science",
            )
        )

for stem in (
    "20181102T041236S859_ovr_scil2_calv2",
    "20181103T041236S523_ovr_scil2_calv2",
):
    for suffix in ("fits", "xml"):
        CORE_FILES.append(
            (
                f"{OVIRS}/data_calibrated/approach/{stem}.{suffix}",
                f"ovirs/sample/{stem}.{suffix}",
                "OVIRS format sample; not the complete fit set",
            )
        )

SHAPE_STEM = "bennu_g_12600mm_alt_obj_0000n00000_v021a"
for suffix in ("bds", "xml"):
    CORE_FILES.append(
        (
            f"{SPICE}/dsk/{SHAPE_STEM}.{suffix}",
            f"shape/{SHAPE_STEM}.{suffix}",
            "public 12.6 m Bennu DSK (v21, close to paper v13 resolution)",
        )
    )

SPICE_FILES = {
    "mk": ("orx_2018_v09.tm", "orx_2018_v09.xml"),
    "lsk": ("naif0012.tls",),
    "pck": ("pck00010.tpc", "bennu_v17.tpc"),
    "fk": ("orx_v14.tf", "orx_shape_v03.tf"),
    "ik": ("orx_otes_v00.ti", "orx_ovirs_v00.ti"),
    "sclk": ("orx_sclkscet_00093.tsc",),
    "spk": (
        "de424.bsp",
        "bennu_refdrmc_v1.bsp",
        "orx_struct_v04.bsp",
        "orx_180801_190302_181218_od077_v1.bsp",
    ),
    "ck": (
        "orx_sc_rel_181029_181104_v02.bc",
        "orx_sc_rel_181105_181111_v02.bc",
    ),
}
for family, names in SPICE_FILES.items():
    for name in names:
        CORE_FILES.append(
            (
                f"{SPICE}/{family}/{name}",
                f"spice/{family}/{name}",
                "SPICE geometry kernel",
            )
        )


def directory_links(url: str) -> list[str]:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        html = response.read().decode("utf-8", errors="replace")
    return [urllib.parse.unquote(x) for x in re.findall(r'href=["\']([^"\']+)', html)]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, relative: str, category: str) -> dict[str, object]:
    destination = DATA / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size:
        print(f"SKIP {relative} ({destination.stat().st_size:,} bytes)")
    else:
        temporary = destination.with_suffix(destination.suffix + ".part")
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        last_error: Exception | None = None
        for attempt in range(1, 4):
            try:
                print(f"GET  {url}")
                with urllib.request.urlopen(request, timeout=120) as response:
                    with temporary.open("wb") as output:
                        while True:
                            block = response.read(1024 * 1024)
                            if not block:
                                break
                            output.write(block)
                temporary.replace(destination)
                last_error = None
                break
            except Exception as exc:  # report URL and retry transient failures
                last_error = exc
                print(f"WARN attempt {attempt}/3 failed: {exc}")
                time.sleep(2 * attempt)
        if last_error is not None:
            raise RuntimeError(f"Failed to download {url}: {last_error}")
    return {
        "category": category,
        "url": url,
        "path": str(destination.relative_to(ROOT)),
        "bytes": destination.stat().st_size,
        "sha256": sha256(destination),
    }


def geometry_files() -> list[tuple[str, str, str]]:
    base = f"{OTES}/geometry/approach/"
    names = directory_links(base)
    pattern = re.compile(r"^2018110[89]T\d{6}S\d{3}_ote_geo\.(?:fits|xml)$")
    selected = sorted({name for name in names if pattern.match(name)})
    if not selected:
        raise RuntimeError("No 2018-11-08/09 OTES geometry products found")
    return [
        (base + name, f"otes/geometry/{name}", "OTES observation geometry")
        for name in selected
    ]


def full_ovirs_files() -> list[tuple[str, str, str]]:
    base = f"{OVIRS}/data_calibrated/approach/"
    names = directory_links(base)
    pattern = re.compile(
        r"^2018110[23]T\d{6}S\d{3}_ovr_scil2_calv2\.(?:fits|xml)$"
    )
    selected = sorted({name for name in names if pattern.match(name)})
    fits_count = sum(name.endswith(".fits") for name in selected)
    if fits_count != 34254:
        raise RuntimeError(
            f"Expected 34,254 OVIRS FITS products but found {fits_count}; "
            "the archive layout may have changed"
        )
    return [
        (base + name, f"ovirs/full/{name}", "complete OVIRS thermal fit set")
        for name in selected
    ]


def write_local_meta_kernel() -> None:
    path = DATA / "spice" / "orex_2018_nov_local.tm"
    kernels = []
    for family, names in SPICE_FILES.items():
        if family == "mk":
            continue
        kernels.extend(f"'{family}/{name}'" for name in names)
    body = (
        "KPL/MK\n\n\\begindata\n\n"
        "PATH_VALUES = ( '.' )\nPATH_SYMBOLS = ( 'ROOT' )\n"
        "KERNELS_TO_LOAD = (\n    $ROOT/"
        + ",\n    $ROOT/".join(kernels)
        + "\n)\n\n\\begintext\n"
        "Local subset for the 2-9 November 2018 Bennu observations.\n"
    )
    path.write_text(body, encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-ovirs",
        action="store_true",
        help="also download the complete two-day OVIRS fit set (~19 GB)",
    )
    parser.add_argument(
        "--confirm-large-download",
        action="store_true",
        help="required safety acknowledgement for --full-ovirs",
    )
    args = parser.parse_args()
    if args.full_ovirs and not args.confirm_large_download:
        parser.error(
            "--full-ovirs needs --confirm-large-download (about 19 GB and "
            "68,508 FITS/XML files)"
        )

    files = list(CORE_FILES)
    files.extend(geometry_files())
    if args.full_ovirs:
        files.extend(full_ovirs_files())

    records = []
    failures = []
    for url, relative, category in files:
        try:
            records.append(download(url, relative, category))
        except Exception as exc:
            print(f"ERROR {exc}")
            failures.append({"url": url, "path": relative, "error": str(exc)})

    write_local_meta_kernel()
    manifest = {
        "description": "Inputs for reproducing the 2019 Bennu thermal analysis",
        "complete_ovirs_requested": args.full_ovirs,
        "files": records,
        "failures": failures,
    }
    manifest_path = DATA / "download_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    total = sum(int(record["bytes"]) for record in records)
    print(f"\nDownloaded/verified {len(records)} files, {total:,} bytes total")
    print(f"Manifest: {manifest_path}")
    if failures:
        raise SystemExit(f"{len(failures)} download(s) failed; see manifest")


if __name__ == "__main__":
    main()
