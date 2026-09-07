"""在 SPCv14 与 SPICE 几何下检查半球坑辐亮度的网格收敛性。"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from pathlib import Path
from time import perf_counter

import numpy as np

from bennu_atpm import ThermoConfig, load_obj, simulate
from reproduce_supplementary_fig2b import (bin_one_degree, geometry_on_model_grid,
                                           process_frames, sampled_mesh)


def main():
    """依次运行 8–200 个坑壁单元，以 200 单元曲线为参考输出差异。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path("output/crater_resolution_convergence.csv"))
    parser.add_argument("--facet-stride", type=int, default=128)
    parser.add_argument("--model-phases", type=int, default=48)
    args = parser.parse_args()
    frame_cache = Path("output/supplementary_fig2b/processed_ovirs_frames.csv")
    rows = process_frames(Path("data/ovirs/full"), frame_cache)
    datasets, _ = bin_one_degree(rows)
    geometry = geometry_on_model_grid(datasets["2018-11-02"], args.model_phases)
    mesh = sampled_mesh(load_obj("data/shape/g_12560mm_spc_obj_0000n00000_v014.obj"),
                        args.facet_stride)
    settings = [(2, 4), (3, 6), (4, 8), (6, 12), (8, 16), (10, 20)]
    results = []
    for theta_bins, azimuth_bins in settings:
        config = ThermoConfig(n_phase=args.model_phases,
                              crater_theta_bins=theta_bins,
                              crater_azimuth_bins=azimuth_bins,
                              facet_chunk_size=16,
                              heliocentric_distance_au=1.041)
        started = perf_counter()
        curve = simulate(mesh, 350.0, config, *geometry).crater_4um
        elapsed = perf_counter()-started
        results.append({"theta_bins": theta_bins, "azimuth_bins": azimuth_bins,
                        "wall_elements": theta_bins*azimuth_bins,
                        "elapsed_s": elapsed, "curve": curve})
        print(f"{theta_bins}x{azimuth_bins}: {elapsed:.2f} s", flush=True)
    # 最高分辨率只作为本组测试的数值参考，并不等于解析真值。
    reference = results[-1]["curve"]
    scale = float(np.mean(reference))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["theta_bins", "azimuth_bins", "wall_elements", "elapsed_s",
                         "mean_crater_4um_W_m-2_um-1_sr-1",
                         "reference_wall_elements",
                         "rms_difference_from_reference_percent",
                         "max_difference_from_reference_percent"])
        for result in results:
            difference = result["curve"]-reference
            writer.writerow([result["theta_bins"], result["azimuth_bins"],
                             result["wall_elements"], result["elapsed_s"],
                             np.mean(result["curve"]),
                             results[-1]["wall_elements"],
                             100*np.sqrt(np.mean(difference**2))/scale,
                             100*np.max(np.abs(difference))/scale])


if __name__ == "__main__":
    main()
