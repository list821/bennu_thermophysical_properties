"""Concurrent, resumable downloader for the OVIRS science FITS used in the paper."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import re
import threading
import time
import urllib.parse
import urllib.request


BASE = "https://sbnarchive.psi.edu/pds4/orex/orex.ovirs/data_calibrated/approach/"
EXPECTED_COUNT = 34254
EXPECTED_SIZE = 558720
USER_AGENT = "Bennu-OVIRS-reproduction/1.0"


def list_products() -> list[str]:
    request = urllib.request.Request(BASE, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        html = response.read().decode("utf-8", errors="replace")
    names = sorted(set(urllib.parse.unquote(name) for name in
                       re.findall(r'href=["\']([^"\']+)', html)
                       if re.match(r"^2018110[23]T\d{6}S\d{3}_ovr_scil2_calv2\.fits$",
                                   urllib.parse.unquote(name))))
    if len(names) != EXPECTED_COUNT:
        raise RuntimeError(f"Expected {EXPECTED_COUNT} FITS products, found {len(names)}")
    return names


def valid_fits(path: Path) -> bool:
    if not path.exists() or path.stat().st_size != EXPECTED_SIZE:
        return False
    with path.open("rb") as stream:
        return stream.read(6) == b"SIMPLE"


def download_one(name: str, destination_dir: Path) -> tuple[str, int, str]:
    destination = destination_dir / name
    if valid_fits(destination):
        return name, destination.stat().st_size, "existing"
    partial = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(1, 6):
        try:
            offset = partial.stat().st_size if partial.exists() else 0
            headers = {"User-Agent": USER_AGENT}
            if offset:
                headers["Range"] = f"bytes={offset}-"
            request = urllib.request.Request(BASE + name, headers=headers)
            with urllib.request.urlopen(request, timeout=180) as response:
                append = offset > 0 and getattr(response, "status", None) == 206
                mode = "ab" if append else "wb"
                with partial.open(mode) as stream:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        stream.write(block)
            if partial.stat().st_size != EXPECTED_SIZE:
                raise IOError(f"size {partial.stat().st_size}, expected {EXPECTED_SIZE}")
            partial.replace(destination)
            if not valid_fits(destination):
                raise IOError("FITS signature/size verification failed")
            return name, destination.stat().st_size, "downloaded"
        except Exception as exc:
            if attempt == 5:
                return name, 0, f"failed: {exc}"
            time.sleep(min(2**attempt, 20))
    return name, 0, "failed"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/ovirs/full"))
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    names = list_products()
    lock = threading.Lock()
    completed = downloaded = existing = failed = total_bytes = 0
    failures = []
    started = time.time()
    print(f"Archive verified: {len(names)} FITS; output={args.output.resolve()}", flush=True)
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(download_one, name, args.output): name for name in names}
        for future in as_completed(futures):
            name, size, status = future.result()
            with lock:
                completed += 1
                total_bytes += size
                if status == "downloaded":
                    downloaded += 1
                elif status == "existing":
                    existing += 1
                else:
                    failed += 1
                    failures.append({"file": name, "error": status})
                if completed % 100 == 0 or completed == len(names):
                    elapsed = max(time.time() - started, 1)
                    rate = total_bytes / elapsed / 1024**2
                    print(f"{completed}/{len(names)} downloaded={downloaded} existing={existing} "
                          f"failed={failed} verified={total_bytes/1024**3:.2f} GiB "
                          f"rate={rate:.1f} MiB/s", flush=True)
    state = {
        "source": BASE,
        "dates": ["2018-11-02", "2018-11-03"],
        "expected_fits": EXPECTED_COUNT,
        "verified_fits": downloaded + existing,
        "bytes": total_bytes,
        "failures": failures,
    }
    (args.output / "download_state.json").write_text(
        json.dumps(state, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    if failures:
        raise SystemExit(f"{len(failures)} files failed; rerun to resume")
    print("OVIRS download complete and verified", flush=True)


if __name__ == "__main__":
    main()
