"""运行 Bennu ATPM 代理正演网格、反演和 100 次 bootstrap。

这里的观测明确是代理数据，并非飞行数据。传入 ``--shape model.obj`` 可替换
内置代理形状；真实 OVIRS 补充图 2b 流程应运行 ``reproduce_supplementary_fig2b.py``。
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from bennu_atpm import (ThermoConfig, crater_fraction_to_rms_deg, generate_grid,
                        load_obj, make_bennu_proxy_mesh, skin_depth_m)


def fit_grid(grid, gammas, fractions, obs4, sigma4, obs14, sigma14):
    """在 Γ–坑覆盖率二维网格上联合计算 4/14 µm χ²。"""
    chi2 = np.empty((len(gammas), len(fractions)))
    for i, gamma in enumerate(gammas):
        model = grid[float(gamma)]
        pred4 = model.smooth_4um[None, :]+fractions[:, None]*(model.crater_4um-model.smooth_4um)
        pred14 = model.smooth_14um[None, :]+fractions[:, None]*(model.crater_14um-model.smooth_14um)
        pred14 /= pred14.mean(axis=1, keepdims=True)
        chi2[i] = (np.sum(((pred4-obs4)/sigma4)**2, axis=1)+
                   np.sum(((pred14-obs14)/sigma14)**2, axis=1))
    index = np.unravel_index(np.argmin(chi2), chi2.shape)
    return chi2, float(gammas[index[0]]), float(fractions[index[1]])


def write_svg(path, phase, obs4, fit4, obs14, fit14):
    """把代理观测和最佳拟合写成双面板 SVG。"""
    width, height, margin = 900, 620, 70
    panels = [(obs4, fit4, "OVIRS 4 um disk radiance proxy"),
              (obs14, fit14, "OTES 14 um relative light curve")]
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#222}.a{stroke:#444}.g{stroke:#ddd}</style>']
    for p, (obs, fit, title) in enumerate(panels):
        top, panel_h = 45+p*285, 225
        left, right, bottom = margin, width-30, top+panel_h
        ymin, ymax = min(obs.min(), fit.min()), max(obs.max(), fit.max())
        pad = max((ymax-ymin)*.1, abs(ymax)*.005, 1e-12)
        ymin, ymax = ymin-pad, ymax+pad
        x = left+phase*(right-left)
        yobs = bottom-(obs-ymin)/(ymax-ymin)*panel_h
        yfit = bottom-(fit-ymin)/(ymax-ymin)*panel_h
        parts.append(f'<text x="{left}" y="{top-14}" font-size="17">{title}</text>')
        for tick in np.linspace(0, 1, 5):
            xx = left+tick*(right-left)
            parts.append(f'<line class="g" x1="{xx}" y1="{top}" x2="{xx}" y2="{bottom}"/>')
            parts.append(f'<text x="{xx}" y="{bottom+20}" text-anchor="middle" font-size="12">{tick:.2f}</text>')
        parts.extend([f'<line class="a" x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}"/>',
                      f'<line class="a" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}"/>'])
        points = " ".join(f"{xx:.2f},{yy:.2f}" for xx, yy in zip(x, yfit))
        parts.append(f'<polyline points="{points}" fill="none" stroke="#d1495b" stroke-width="2.2"/>')
        parts.extend(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="2.8" fill="#1769aa"/>'
                     for xx, yy in zip(x[::4], yobs[::4]))
    parts.append(f'<text x="{width/2}" y="{height-8}" text-anchor="middle">Rotational phase</text></svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def main():
    """生成代理观测、执行网格反演和 bootstrap，并写出全部结果。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shape", type=Path, help="optional triangular/polygonal Bennu OBJ")
    parser.add_argument("--output", type=Path, default=Path("output/python_bennu"))
    parser.add_argument("--bootstrap", type=int, default=100)
    parser.add_argument("--quick", action="store_true", help="50-unit Gamma spacing and 64 phases")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    config = ThermoConfig(n_phase=64 if args.quick else 96)
    mesh = load_obj(args.shape) if args.shape else make_bennu_proxy_mesh()
    gammas = np.arange(0., 601., 50. if args.quick else 10.)
    # 0.01 网格精确包含论文报告的坑覆盖率 0.77。
    fractions = np.round(np.arange(0., 1.0001, .01), 2)
    print(f"Mesh: {mesh.source}; {len(mesh.faces)} facets")
    print(f"Calculating {len(gammas)} thermal-inertia models ...")
    grid = generate_grid(mesh, gammas, config)

    truth_gamma, truth_fraction = (350., .77)
    truth = grid[truth_gamma]
    obs4, obs14 = truth.mixed(truth_fraction)
    obs14 /= obs14.mean()
    # 这是人为代理误差，不是飞行数据不确定度；噪声大小只用于复现论文 bootstrap
    # 宽度的量级（Γ 约 20、坑覆盖率约 0.04）。
    sigma4 = np.full_like(obs4, .010*obs4.mean())
    sigma14 = np.full_like(obs14, .0035)
    chi2, best_gamma, best_fraction = fit_grid(grid, gammas, fractions, obs4, sigma4, obs14, sigma14)

    rng = np.random.default_rng(2019)
    bootstrap = []
    for _ in range(args.bootstrap):
        sample4 = obs4+rng.normal(0, sigma4)
        sample14 = obs14+rng.normal(0, sigma14)
        _, gamma, fraction = fit_grid(grid, gammas, fractions, sample4, sigma4, sample14, sigma14)
        bootstrap.append((gamma, fraction, crater_fraction_to_rms_deg(fraction)))
    boot = np.asarray(bootstrap)
    ddof = 1 if len(boot) > 1 else 0
    sigmas = np.std(boot, axis=0, ddof=ddof) if len(boot) else np.full(3, np.nan)
    best4, best14 = grid[best_gamma].mixed(best_fraction)
    best14 /= best14.mean()

    with (args.output/"best_fit_curves.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["phase", "proxy_ovirs_4um", "fit_ovirs_4um",
                         "proxy_otes_14um_relative", "fit_otes_14um_relative"])
        writer.writerows(zip(truth.phase, obs4, best4, obs14, best14))
    with (args.output/"chi2_grid.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["gamma", "crater_fraction", "rms_slope_deg", "chi2"])
        for i, gamma in enumerate(gammas):
            for j, fraction in enumerate(fractions):
                writer.writerow([gamma, fraction, crater_fraction_to_rms_deg(fraction), chi2[i, j]])
    with (args.output/"bootstrap.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["thermal_inertia", "crater_fraction", "rms_slope_deg"])
        writer.writerows(bootstrap)
    write_svg(args.output/"lightcurve_comparison.svg", truth.phase, obs4, best4, obs14, best14)

    summary = {
        "status": "publication-constrained proxy reproduction; not archived flight-data fitting",
        "shape": mesh.source, "facet_count": int(len(mesh.faces)),
        "paper_targets": {"thermal_inertia": [350., 20.], "rms_slope_deg": [43., 1.],
                          "crater_fraction": [.77, .04]},
        "recovered": {"thermal_inertia": best_gamma,
                      "bootstrap_sigma_thermal_inertia": float(sigmas[0]),
                      "crater_fraction": best_fraction,
                      "bootstrap_sigma_crater_fraction": float(sigmas[1]),
                      "rms_slope_deg": crater_fraction_to_rms_deg(best_fraction),
                      "bootstrap_sigma_rms_slope_deg": float(sigmas[2])},
        "physics_checks": {"nominal_skin_depth_cm": 100*skin_depth_m(350, config),
                           "largest_boundary_residual_W_m2": float(max(x.max_boundary_residual_w_m2 for x in grid.values()))},
        "model": {"conduction": "periodic semi-infinite 1-D Fourier solution",
                  "boundary": "absorbed sunlight = epsilon*sigma*T^4 + downward conduction",
                  "roughness": "fractional energy-conserving hemispherical-crater microfacets",
                  "included": ["crater aperture visibility", "crater self-heating",
                               "multiple-scattered sunlight", "Gaussian OVIRS sample response"],
                  "omissions": ["global shape shadow ray tracing", "global shape mutual heating",
                                "SPICE geometry", "full measured instrument response"]},
        "grid": {"gamma_min": float(gammas.min()), "gamma_max": float(gammas.max()),
                 "gamma_step": float(gammas[1]-gammas[0]), "fraction_step": .01,
                 "bootstrap_trials": args.bootstrap}}
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
