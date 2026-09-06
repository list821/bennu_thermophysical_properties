"""Fit the public calibrated OSIRIS-REx thermal observations with the Python model.

This is a real-data inversion, but not yet a strict reproduction of the 2019
paper: the public v21 DSK is not read without SPICE/DSK support, so the built-in
proxy shape is used unless --shape supplies an OBJ.  Each observing day's
rotation zero point is fitted as a nuisance parameter.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from bennu_atpm import (ThermoConfig, crater_fraction_to_rms_deg, generate_grid,
                        load_obj, make_bennu_proxy_mesh)
from pds_observations import (BinnedLightcurve, load_otes_lightcurves,
                              load_ovirs_lightcurves, load_ovirs_samples)


def shifted_chi2(prediction: np.ndarray, observation: BinnedLightcurve) -> tuple[float, int]:
    valid = np.isfinite(observation.value) & np.isfinite(observation.sigma) & (observation.sigma > 0)
    if np.count_nonzero(valid) < 4:
        return np.inf, 0
    best, best_shift = np.inf, 0
    for shift in range(len(prediction)):
        candidate = np.roll(prediction, shift)
        value = float(np.sum(((candidate[valid] - observation.value[valid]) /
                              observation.sigma[valid]) ** 2))
        if value < best:
            best, best_shift = value, shift
    return best, best_shift


def fit_grid(grids: dict[str, dict[float, object]], observations: list[BinnedLightcurve],
             gammas: np.ndarray, fractions: np.ndarray) -> tuple[np.ndarray, float, float]:
    chi2 = np.zeros((len(gammas), len(fractions)), float)
    for i, gamma in enumerate(gammas):
        for j, fraction in enumerate(fractions):
            total = 0.0
            for observation in observations:
                model = grids[observation.instrument][float(gamma)]
                curve4, curve14 = model.mixed(float(fraction))
                prediction = curve4 if observation.instrument == "OVIRS" else curve14
                prediction = prediction * model.projected_area_m2
                prediction = prediction / np.mean(prediction)
                value, _ = shifted_chi2(prediction, observation)
                total += value
            chi2[i, j] = total
    location = np.unravel_index(np.nanargmin(chi2), chi2.shape)
    return chi2, float(gammas[location[0]]), float(fractions[location[1]])


def best_details(grids, observations, gamma, fraction):
    output = []
    for observation in observations:
        model = grids[observation.instrument][gamma]
        curve4, curve14 = model.mixed(fraction)
        prediction = curve4 if observation.instrument == "OVIRS" else curve14
        prediction = prediction * model.projected_area_m2
        prediction = prediction / np.mean(prediction)
        value, shift = shifted_chi2(prediction, observation)
        output.append({
            "instrument": observation.instrument,
            "source": observation.source,
            "chi2": value,
            "valid_bins": int(np.count_nonzero(np.isfinite(observation.value))),
            "phase_shift_bins": shift,
            "phase_shift_cycles": shift / len(prediction),
            "model": np.roll(prediction, shift),
        })
    return output


def write_curve_csv(path: Path, observations, details):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["instrument", "source", "phase", "observed_relative_radiance",
                         "sigma", "sample_count", "best_model_relative_radiance"])
        for observation, detail in zip(observations, details):
            for row in zip(observation.phase, observation.value, observation.sigma,
                           observation.count, detail["model"]):
                writer.writerow([observation.instrument, observation.source, *row])


def write_svg(path: Path, observations, details):
    width, panel_h, margin = 980, 235, 65
    height = 55 + len(observations) * (panel_h + 50)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#222}.grid{stroke:#ddd}.axis{stroke:#333}</style>']
    for panel, (obs, detail) in enumerate(zip(observations, details)):
        top = 45 + panel * (panel_h + 50)
        bottom, left, right = top + panel_h, margin, width - 25
        valid = np.isfinite(obs.value)
        all_values = np.r_[obs.value[valid], detail["model"][valid]]
        ymin, ymax = float(np.min(all_values)), float(np.max(all_values))
        pad = max((ymax - ymin) * .12, .005)
        ymin, ymax = ymin - pad, ymax + pad
        parts.append(f'<text x="{left}" y="{top-13}" font-size="16">{obs.instrument}: {obs.source}</text>')
        for tick in np.linspace(0, 1, 5):
            x = left + tick * (right - left)
            parts.append(f'<line class="grid" x1="{x}" y1="{top}" x2="{x}" y2="{bottom}"/>')
            parts.append(f'<text x="{x}" y="{bottom+19}" text-anchor="middle" font-size="12">{tick:.2f}</text>')
        parts.extend([f'<line class="axis" x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}"/>',
                      f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}"/>'])
        x = left + obs.phase * (right - left)
        ymodel = bottom - (detail["model"] - ymin) / (ymax - ymin) * panel_h
        points = " ".join(f"{xx:.2f},{yy:.2f}" for xx, yy in zip(x, ymodel))
        parts.append(f'<polyline points="{points}" fill="none" stroke="#d1495b" stroke-width="2"/>')
        yobs = bottom - (obs.value - ymin) / (ymax - ymin) * panel_h
        parts.extend(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="2.2" fill="#1769aa"/>'
                     for xx, yy, ok in zip(x, yobs, valid) if ok)
    parts.append(f'<text x="{width/2}" y="{height-8}" text-anchor="middle">Rotational phase (daily zero point fitted)</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data"))
    parser.add_argument("--shape", type=Path, help="optional Bennu OBJ; otherwise proxy shape")
    parser.add_argument("--output", type=Path, default=Path("output/real_bennu"))
    parser.add_argument("--gamma-step", type=float, default=10.0)
    parser.add_argument("--fraction-step", type=float, default=0.01)
    parser.add_argument("--phase-bins", type=int, default=96)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    otes, otes_diagnostics = load_otes_lightcurves(
        args.data / "otes/science", args.data / "otes/geometry", n_bins=args.phase_bins
    )
    ovirs_directory = args.data / "ovirs/full"
    if not ovirs_directory.exists():
        ovirs_directory = args.data / "ovirs/sample"
    ovirs, ovirs_diagnostics = load_ovirs_lightcurves(ovirs_directory, n_bins=args.phase_bins)
    ovirs_samples = load_ovirs_samples(ovirs_directory)
    observations = otes + ovirs
    if not observations:
        raise SystemExit("No usable real thermal light curves were found")

    mesh = load_obj(args.shape) if args.shape else make_bennu_proxy_mesh()
    gammas = np.arange(0.0, 600.0 + args.gamma_step / 2, args.gamma_step)
    fractions = np.arange(0.0, 1.0 + args.fraction_step / 2, args.fraction_step)
    base = ThermoConfig(n_phase=args.phase_bins)
    configs = {
        "OTES": replace(base, heliocentric_distance_au=1.0195, phase_angle_deg=5.1),
        "OVIRS": replace(base, heliocentric_distance_au=1.039, phase_angle_deg=4.8),
    }
    grids = {}
    for instrument in sorted({item.instrument for item in observations}):
        print(f"Computing {instrument} model grid: {len(gammas)} thermal inertias")
        grids[instrument] = generate_grid(mesh, gammas, configs[instrument])
    chi2, gamma, fraction = fit_grid(grids, observations, gammas, fractions)
    details = best_details(grids, observations, gamma, fraction)

    with (args.output / "real_chi2_grid.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["thermal_inertia", "crater_fraction", "rms_slope_deg", "chi2"])
        for i, g in enumerate(gammas):
            for j, f in enumerate(fractions):
                writer.writerow([g, f, crater_fraction_to_rms_deg(f), chi2[i, j]])
    write_curve_csv(args.output / "real_lightcurves.csv", observations, details)
    write_svg(args.output / "real_lightcurve_fit.svg", observations, details)

    min_chi2 = float(np.min(chi2))
    dof = sum(int(np.count_nonzero(np.isfinite(item.value))) for item in observations) - 2 - len(observations)
    summary = {
        "status": "real calibrated PDS observations fitted with simplified Python ATPM",
        "result_is_strict_2019_reproduction": False,
        "best_fit": {
            "thermal_inertia_J_m-2_K-1_s-1/2": gamma,
            "crater_fraction": fraction,
            "rms_slope_deg": crater_fraction_to_rms_deg(fraction),
            "chi2": min_chi2,
            "degrees_of_freedom_approx": dof,
            "reduced_chi2_approx": min_chi2 / max(dof, 1),
        },
        "data": {
            "OTES": otes_diagnostics,
            "OVIRS_lightcurve": ovirs_diagnostics,
            "OVIRS_extracted_frames": ovirs_samples,
            "OVIRS_used_in_fit": bool(ovirs),
        },
        "model": {
            "shape": mesh.source,
            "facet_count": int(len(mesh.faces)),
            "gamma_step": args.gamma_step,
            "fraction_step": args.fraction_step,
            "phase_bins": args.phase_bins,
            "nuisance_parameters": "independent cyclic phase shift for each observing day",
            "limitations": [
                "proxy shape unless --shape OBJ is supplied",
                "mean observation geometry rather than per-facet SPICE geometry",
                "no global ray-traced shadowing or iterative self-heating",
                "monochromatic band averages rather than full instrument response",
            ],
        },
        "per_curve_fit": [{k: v for k, v in item.items() if k != "model"} for item in details],
    }
    (args.output / "real_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
