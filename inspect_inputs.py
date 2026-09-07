"""以人可读方式检查当前 ATPM 所需的 OVIRS、形状、太阳谱和 SPICE 输入。"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

from bennu_atpm import load_obj
from download_spice_kernels import KERNELS
from pds_observations import extract_ovirs_thermal_band, read_ovirs_product


def main() -> None:
    """读取一个实际样本并打印数组形状、单位、元数据和运行所需文件状态。"""
    # Windows 中文系统的默认 GBK 不能表示 µ 等字符；统一为 UTF-8 便于 VS Code 显示。
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ovirs", type=Path,
                        default=Path("data/ovirs/sample/20181102T041236S859_ovr_scil2_calv2.fits"))
    parser.add_argument("--shape", type=Path,
                        default=Path("data/shape/g_12560mm_spc_obj_0000n00000_v014.obj"))
    parser.add_argument("--kernels", type=Path,
                        default=Path(r"D:\ATPM\data\spice\kernels"))
    parser.add_argument("--processed", type=Path,
                        default=Path("output/supplementary_fig2b/processed_ovirs_frames.csv"))
    args = parser.parse_args()

    product = read_ovirs_product(args.ovirs)
    thermal, sigma, samples = extract_ovirs_thermal_band(product)
    print("OVIRS FITS")
    print(f"  文件: {args.ovirs.resolve()}")
    print(f"  观测时间: {product['mid_obs'].isoformat()}")
    print(f"  radiance: {product['radiance'].shape}, W cm^-2 um^-1 sr^-1")
    print(f"  wavelength_um: {product['wavelength_um'].shape}, um")
    print(f"  uncertainty: {product['uncertainty'].shape}, 与 radiance 同单位")
    print(f"  quality: {product['quality'].shape}, 位掩码")
    print(f"  PDS FILL_FAC: {float(product['fov_fill_factor']):.6g}（仅诊断，不再相除）")
    print(f"  扣反射后的 4 um 热辐亮度: {thermal:.9g} ± {sigma:.3g}; 样本数={samples}")

    solar_path = Path("data/ovirs/orexsolarflux.csv")
    solar = np.loadtxt(solar_path, delimiter=",")
    print("\n太阳光谱 CSV")
    print(f"  文件: {solar_path.resolve()}")
    print(f"  形状: {solar.shape}; 第1列=波长，第2列=光谱能量密度，第3列=不确定度，第4列=平滑分辨率")
    print(f"  波长范围: {solar[:, 0].min():.6g} 至 {solar[:, 0].max():.6g} um")

    mesh = load_obj(args.shape)
    print("\nSPCv14 OBJ")
    print(f"  文件: {args.shape.resolve()}")
    print(f"  顶点数: {len(mesh.vertices)}; 三角面元数: {len(mesh.faces)}")
    print(f"  面元面积范围: {mesh.areas.min():.6g} 至 {mesh.areas.max():.6g} m^2")
    print("  OBJ 的 v 行是三维顶点；f 行是从 1 开始的顶点索引。PDS 坐标为 km，读取时转成 m。")

    missing = [relative for relative in KERNELS if not (args.kernels / relative).exists()]
    print("\nSPICE kernels")
    print(f"  根目录: {args.kernels}")
    print(f"  需要 {len(KERNELS)} 个，缺失 {len(missing)} 个")
    for relative in missing:
        print(f"  缺失: {relative}")

    print("\n逐帧处理缓存")
    if args.processed.exists():
        with args.processed.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            first = next(reader, None)
            print(f"  文件: {args.processed.resolve()}")
            print(f"  字段数: {len(reader.fieldnames or [])}")
            print("  字段: " + ", ".join(reader.fieldnames or []))
            if first:
                print(f"  首帧: {first.get('utc')}; thermal_radiance={first.get('thermal_radiance')}")
    else:
        print("  尚未生成；第一次运行主程序后会出现。")


if __name__ == "__main__":
    main()
