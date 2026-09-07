"""把 NAIF type-2 DSK 板片形状转换为普通三角形 OBJ。"""

from __future__ import annotations

import argparse
from pathlib import Path

import spiceypy as spice


def main():
    """读取 DSK 首个段的顶点和板片，换算单位后写成 OBJ。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    handle = spice.dasopr(str(args.input.resolve()))
    try:
        descriptor = spice.dlabfs(handle)
        vertex_count, plate_count = spice.dskz02(handle, descriptor)
        vertices = spice.dskv02(handle, descriptor, 1, vertex_count)
        plates = spice.dskp02(handle, descriptor, 1, plate_count)
    finally:
        spice.dascls(handle)
    # DSK 顶点单位为 km；转换为 m 以匹配 bennu_atpm.py。
    with args.output.open("w", encoding="ascii", newline="\n") as stream:
        stream.write(f"# Converted from {args.input.name}; {vertex_count} vertices, {plate_count} plates\n")
        for x, y, z in vertices:
            stream.write(f"v {x*1000:.9f} {y*1000:.9f} {z*1000:.9f}\n")
        for a, b, c in plates:  # DSK 与 OBJ 的面索引在这里都从 1 开始。
            stream.write(f"f {int(a)} {int(b)} {int(c)}\n")
    print(f"Wrote {args.output}: {vertex_count} vertices, {plate_count} triangular facets")


if __name__ == "__main__":
    main()
