"""下载计算 2018 年 11 月观测几何所需的最小官方 SPICE 核集合。"""
from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import urlopen


BASE = "https://naif.jpl.nasa.gov/pub/naif/pds/pds4/orex/orex_spice/spice_kernels"
KERNELS = (
    "lsk/naif0012.tls",
    "pck/pck00010.tpc",
    # 必须匹配 SPCv14 OBJ 的经度约定；加载 bennu_v17 会让 v14 表面绕自转轴偏转。
    "pck/bennu_v14.tpc",
    "fk/orx_v14.tf",
    "fk/orx_shape_v03.tf",
    "ik/orx_ovirs_v00.ti",
    "sclk/orx_sclkscet_00093.tsc",
    "spk/de424.bsp",
    "spk/bennu_refdrmc_v1.bsp",
    "spk/orx_struct_v04.bsp",
    "spk/orx_180801_190302_181218_od077_v1.bsp",
    # 覆盖两段观测的重建航天器姿态 CK。
    "ck/orx_sc_rel_181029_181104_v02.bc",
)


def download(destination: Path) -> None:
    """跳过已存在核文件，其余逐个下载到保持 NAIF 目录结构的目标位置。"""
    for relative in KERNELS:
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and target.stat().st_size > 0:
            print(f"exists {relative} ({target.stat().st_size:,} bytes)", flush=True)
            continue
        temporary = target.with_suffix(target.suffix + ".part")
        print(f"download {relative}", flush=True)
        with urlopen(f"{BASE}/{relative}", timeout=120) as response, temporary.open("wb") as out:
            total = int(response.headers.get("Content-Length", 0))
            copied = 0
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                out.write(block)
                copied += len(block)
                if copied % (25 * 1024 * 1024) < len(block):
                    print(f"  {copied/1024**2:.0f}/{total/1024**2:.0f} MiB", flush=True)
        temporary.replace(target)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/spice/kernels"))
    download(parser.parse_args().output)
