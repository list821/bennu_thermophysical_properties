"""用公开 OVIRS 数据重建 4 µm 光变和论文补充图 2b。

输入为 PDS 公开的 L2 v2 标定 FITS。程序在约 2.1 µm 处缩放太阳光谱形状，
从 4 µm 测量中减去反射分量，再按 Bennu 自转相位每 1° 分箱，以近似作者未
公开的 L3a 序列。反演固定 RMS 粗糙度 43°（半球坑覆盖率 0.77），扫描热惯量。
"""

from __future__ import annotations

import argparse
import csv
import gc
from dataclasses import replace
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np

from bennu_atpm import Mesh, ThermoConfig, load_obj, make_bennu_proxy_mesh, simulate
from pds_observations import extract_ovirs_thermal_band, read_ovirs_product
from spice_geometry import KERNEL_SET_ID, frame_geometry


DATE_CONFIG = {
    "2018-11-02": {"heliocentric_distance_au": 1.041, "phase_angle_deg": 5.1},
    "2018-11-03": {"heliocentric_distance_au": 1.037, "phase_angle_deg": 4.5},
}


def process_frames(input_dir: Path, cache: Path) -> list[dict[str, object]]:
    """读取两天的 OVIRS FITS、扣除反射光并缓存逐帧热辐亮度。"""
    fits_files = sorted(input_dir.glob("2018110[23]T*_ovr_scil2_calv2.fits"))
    if not fits_files:
        raise FileNotFoundError(f"No OVIRS science FITS in {input_dir}")
    if cache.exists():
        with cache.open("r", newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        # 论文从辐射定标光谱开始，并未要求再除以 FILL_FAC。这里只把 FILL_FAC
        # 留作诊断；反演使用检测器口径辐亮度与孔径积分模型直接比较。
        cache_is_current = (bool(rows) and
                            rows[0].get("band_definition") == "3.98-4.02 um" and
                            rows[0].get("solar_spectrum") == "PDS orexsolarflux.csv")
        if cache_is_current and len(rows) >= int(len(fits_files) * 0.8):
            print(f"Using cached processed frames: {len(rows)}", flush=True)
            for row in rows:
                for key in ("unix_s", "thermal_radiance", "sigma", "fill_factor"):
                    row[key] = float(row[key])
                raw = float(row.get("thermal_radiance_detector", row["thermal_radiance"]))
                raw_sigma = float(row.get("sigma_detector", row["sigma"]))
                fill = float(row["fill_factor"])
                row["thermal_radiance_detector"] = raw
                row["sigma_detector"] = raw_sigma
                row["thermal_radiance"] = raw
                row["sigma"] = raw_sigma
                row["disk_equivalent_radiance_diagnostic"] = raw/fill
                row["radiance_correction"] = "none_PDS_FILL_FAC_is_diagnostic_only"
            # 把正确的数据口径写回缓存，避免脱离本程序查看 CSV 时产生误解。
            with cache.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            return rows
    rows = []
    print(f"Processing {len(fits_files)} calibrated OVIRS FITS", flush=True)
    for index, path in enumerate(fits_files, 1):
        product = read_ovirs_product(path)
        thermal, sigma, spectral_samples = extract_ovirs_thermal_band(product)
        fill = float(product["fov_fill_factor"])
        if np.isfinite(thermal) and thermal > 0 and np.isfinite(fill) and fill > 0:
            timestamp = product["mid_obs"]
            rows.append({
                "file": path.name,
                "date": timestamp.date().isoformat(),
                "utc": timestamp.isoformat(),
                "unix_s": timestamp.timestamp(),
                "thermal_radiance_detector": thermal,
                "sigma_detector": max(sigma, thermal * 0.002),
                "thermal_radiance": thermal,
                "sigma": max(sigma, thermal * 0.002),
                "disk_equivalent_radiance_diagnostic": thermal/fill,
                "fill_factor": fill,
                "spectral_samples": spectral_samples,
                "radiance_correction": "none_PDS_FILL_FAC_is_diagnostic_only",
                "band_definition": "3.98-4.02 um",
                "solar_spectrum": "PDS orexsolarflux.csv",
            })
        if index % 1000 == 0 or index == len(fits_files):
            print(f"Processed {index}/{len(fits_files)}; accepted={len(rows)}", flush=True)
    cache.parent.mkdir(parents=True, exist_ok=True)
    with cache.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


SPICE_FLOAT_FIELDS = (
    "et_s", "sun_x", "sun_y", "sun_z", "observer_x", "observer_y", "observer_z",
    "heliocentric_distance_au", "spacecraft_distance_km", "phase_angle_deg",
    "subsolar_latitude_deg",
    "boresight_x", "boresight_y", "boresight_z", "fov_half_angle_rad",
)


def attach_spice_geometry(rows: list[dict[str, object]], kernel_root: Path,
                          cache: Path) -> list[dict[str, object]]:
    """按每帧 UTC 计算并附加太阳、航天器、视轴与 FOV 的 SPICE 几何。"""
    if rows and all(row.get("spice_kernel_set") == KERNEL_SET_ID for row in rows):
        for row in rows:
            for key in SPICE_FLOAT_FIELDS:
                row[key] = float(row[key])
        print(f"Using cached per-frame SPICE geometry: {len(rows)}", flush=True)
        return rows
    print(f"Evaluating SPICE geometry for {len(rows)} OVIRS frames", flush=True)
    geometry = frame_geometry((str(row["utc"]) for row in rows), kernel_root)
    for index, row in enumerate(rows):
        row.update({
            "spice_frame": "IAU_BENNU",
            "spice_aberration": "LT+S",
            "spice_kernel_set": KERNEL_SET_ID,
            "et_s": geometry["et_s"][index],
            "sun_x": geometry["sun_direction"][index, 0],
            "sun_y": geometry["sun_direction"][index, 1],
            "sun_z": geometry["sun_direction"][index, 2],
            "observer_x": geometry["observer_direction"][index, 0],
            "observer_y": geometry["observer_direction"][index, 1],
            "observer_z": geometry["observer_direction"][index, 2],
            "heliocentric_distance_au": geometry["heliocentric_distance_au"][index],
            "spacecraft_distance_km": geometry["spacecraft_distance_km"][index],
            "phase_angle_deg": geometry["phase_angle_deg"][index],
            "subsolar_latitude_deg": geometry["subsolar_latitude_deg"][index],
            "boresight_x": geometry["ovirs_boresight_direction"][index, 0],
            "boresight_y": geometry["ovirs_boresight_direction"][index, 1],
            "boresight_z": geometry["ovirs_boresight_direction"][index, 2],
            "fov_half_angle_rad": geometry["ovirs_fov_half_angle_rad"][index],
        })
    with cache.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def bin_one_degree(rows: list[dict[str, object]], period_hours: float = 4.296061):
    """按 4.296061 h 自转周期将逐帧量合并到 360 个一度相位箱。"""
    epoch = min(float(row["unix_s"]) for row in rows)
    period_s = period_hours * 3600.0
    datasets = {}
    for day in DATE_CONFIG:
        selected = [row for row in rows if row["date"] == day]
        phase = np.mod((np.asarray([float(row["unix_s"]) for row in selected]) - epoch) / period_s, 1.0)
        value = np.asarray([float(row["thermal_radiance"]) for row in selected])
        detector_value = np.asarray([float(row["thermal_radiance_detector"])
                                     for row in selected])
        disk_equivalent_value = np.asarray(
            [float(row.get("disk_equivalent_radiance_diagnostic",
                           float(row["thermal_radiance_detector"])/float(row["fill_factor"])))
             for row in selected])
        frame_sigma = np.asarray([float(row["sigma"]) for row in selected])
        bin_index = np.floor(phase * 360).astype(int)
        centers = (np.arange(360) + 0.5) / 360.0
        means = np.full(360, np.nan)
        errors = np.full(360, np.nan)
        counts = np.zeros(360, int)
        detector_means = np.full(360, np.nan)
        disk_equivalent_means = np.full(360, np.nan)
        fill_means = np.full(360, np.nan)
        sun_direction = np.full((360, 3), np.nan)
        observer_direction = np.full((360, 3), np.nan)
        heliocentric_distance = np.full(360, np.nan)
        observer_distance = np.full(360, np.nan)
        boresight_direction = np.full((360, 3), np.nan)
        fov_half_angle = np.full(360, np.nan)
        for index in range(360):
            mask = bin_index == index
            samples = value[mask]
            sample_errors = frame_sigma[mask]
            if not len(samples):
                continue
            median = np.median(samples)
            mad = 1.4826 * np.median(np.abs(samples - median))
            keep = np.ones(len(samples), bool) if mad == 0 else np.abs(samples - median) <= 5 * mad
            samples = samples[keep]
            sample_errors = sample_errors[keep]
            counts[index] = len(samples)
            means[index] = float(np.mean(samples))
            scatter_sem = (float(np.std(samples, ddof=1) / np.sqrt(len(samples)))
                           if len(samples) > 1 else 0.0)
            spectral_sem = float(np.sqrt(np.sum(sample_errors**2)) / max(len(sample_errors), 1))
            # 0.5% 误差下限防止大量光谱把尚未建模的形状细节变成不现实的小误差。
            errors[index] = max(scatter_sem, spectral_sem, means[index] * 0.005)
            kept_indices = np.flatnonzero(mask)[keep]
            detector_means[index] = np.mean(detector_value[kept_indices])
            disk_equivalent_means[index] = np.mean(disk_equivalent_value[kept_indices])
            fill_means[index] = np.mean([float(selected[j]["fill_factor"])
                                         for j in kept_indices])
            if selected and "sun_x" in selected[0]:
                sun_mean = np.mean([[float(selected[j][f"sun_{axis}"]) for axis in "xyz"]
                                    for j in kept_indices], axis=0)
                observer_mean = np.mean(
                    [[float(selected[j][f"observer_{axis}"]) for axis in "xyz"]
                     for j in kept_indices], axis=0)
                sun_direction[index] = sun_mean/np.linalg.norm(sun_mean)
                observer_direction[index] = observer_mean/np.linalg.norm(observer_mean)
                heliocentric_distance[index] = np.mean(
                    [float(selected[j]["heliocentric_distance_au"]) for j in kept_indices])
                observer_distance[index] = np.mean(
                    [float(selected[j]["spacecraft_distance_km"]) for j in kept_indices])
                if "boresight_x" in selected[0]:
                    bore_mean = np.mean(
                        [[float(selected[j][f"boresight_{axis}"]) for axis in "xyz"]
                         for j in kept_indices], axis=0)
                    boresight_direction[index] = bore_mean/np.linalg.norm(bore_mean)
                    fov_half_angle[index] = np.mean(
                        [float(selected[j]["fov_half_angle_rad"]) for j in kept_indices])
        datasets[day] = {"phase": centers, "value": means, "sigma": errors, "count": counts,
                         "detector_value": detector_means,
                         "disk_equivalent_diagnostic": disk_equivalent_means,
                         "pds_fill_factor": fill_means,
                         "sun_direction": sun_direction,
                         "observer_direction": observer_direction,
                         "heliocentric_distance_au": heliocentric_distance,
                         "observer_distance_km": observer_distance,
                         "boresight_direction": boresight_direction,
                         "fov_half_angle_rad": fov_half_angle}
    return datasets, epoch


def cyclic_interpolate(curve: np.ndarray, phase: np.ndarray) -> np.ndarray:
    """在周期边界连续地把模型曲线插值到观测相位。"""
    model_phase = np.arange(len(curve)) / len(curve)
    return np.interp(np.mod(phase, 1.0), np.r_[model_phase, 1.0], np.r_[curve, curve[0]])


def geometry_on_model_grid(data, n_phase: int):
    """把一度分箱后的逐帧 SPICE 几何周期插值到热模型相位网格。"""
    target = np.arange(n_phase)/n_phase
    valid = np.isfinite(data["heliocentric_distance_au"])
    source_phase = data["phase"][valid]
    if np.count_nonzero(valid) < 3:
        return None

    def interpolate(values):
        values = np.asarray(values)[valid]
        xp = np.r_[source_phase-1, source_phase, source_phase+1]
        if values.ndim == 1:
            return np.interp(target, xp, np.tile(values, 3))
        return np.column_stack([np.interp(target, xp, np.tile(values[:, i], 3))
                                for i in range(values.shape[1])])

    sun = interpolate(data["sun_direction"])
    observer = interpolate(data["observer_direction"])
    sun /= np.linalg.norm(sun, axis=1)[:, None]
    observer /= np.linalg.norm(observer, axis=1)[:, None]
    boresight = interpolate(data["boresight_direction"])
    boresight /= np.linalg.norm(boresight, axis=1)[:, None]
    return (sun, observer, interpolate(data["heliocentric_distance_au"]),
            interpolate(data["observer_distance_km"]), boresight,
            interpolate(data["fov_half_angle_rad"]))


def sampled_mesh(mesh: Mesh, stride: int) -> Mesh:
    """为热惯量扫描均匀抽取形状面元，并按步长补偿代表面积。"""
    stride = max(1, int(stride))
    index = np.arange(0, len(mesh.faces), stride)
    # 每个保留面元代表邻近被跳过的面元。面积比例在盘平均辐亮度中会约去，
    # 但对 FOV 立体角和扫描网格上的近似全局辐射交换不可缺少。
    return Mesh(mesh.vertices, mesh.faces[index], mesh.centers[index], mesh.normals[index],
                mesh.areas[index]*stride, f"{mesh.source}; uniform plate stride {stride}")


def fit_gamma(mesh, datasets, gammas, roughness_fraction, model_phases,
              crater_theta_bins=8, crater_azimuth_bins=16,
              facet_chunk_size=64, self_heating_iterations=4,
              global_shadowing=True, global_self_heating=True,
              view_factor_cache: Path | None = None):
    """逐个 Γ 正演并以两个观测日的加权残差联合搜索最小 χ²。"""
    fit_cache = None
    if view_factor_cache is not None:
        fit_cache = str(view_factor_cache.with_name(
            f"{view_factor_cache.stem}_{len(mesh.faces)}facets{view_factor_cache.suffix}"))
    config0 = ThermoConfig(n_phase=model_phases, crater_theta_bins=crater_theta_bins,
                           crater_azimuth_bins=crater_azimuth_bins,
                           facet_chunk_size=facet_chunk_size,
                           crater_self_heating_iterations=self_heating_iterations,
                           global_shadowing=global_shadowing,
                           global_self_heating=global_self_heating,
                           global_view_factor_cache=fit_cache)
    model = {day: {} for day in datasets}
    for day in datasets:
        config = replace(config0, **DATE_CONFIG[day])
        geometry = geometry_on_model_grid(datasets[day], model_phases)
        print(f"Computing {day} models ({len(gammas)} Gamma values)", flush=True)
        for gamma in gammas:
            if geometry is None:
                result = simulate(mesh, float(gamma), config)
            else:
                result = simulate(mesh, float(gamma), config, *geometry)
            curve4 = result.mixed_detector(roughness_fraction)
            # 模型核心输出 W m⁻²；论文图纵轴使用 W cm⁻²，故除以 10⁴。
            model[day][float(gamma)] = curve4 / 1.0e4
            gc.collect()
    # SPICE 与 SPC 形状共享天体固定经度约定。使用 SPICE 时再拟合任意相位平移
    # 会丢弃核文件提供的信息；只有缺少 SPICE 几何时才保留旧的相位扫描。
    using_spice = all(np.any(np.isfinite(data["heliocentric_distance_au"]))
                      for data in datasets.values())
    shifts = np.array([0.0]) if using_spice else np.arange(360)/360.0
    chi2 = np.zeros((len(gammas), len(shifts)))
    for gi, gamma in enumerate(gammas):
        for si, shift in enumerate(shifts):
            total = 0.0
            for day, observed in datasets.items():
                valid = np.isfinite(observed["value"]) & np.isfinite(observed["sigma"])
                prediction = cyclic_interpolate(model[day][float(gamma)], observed["phase"] - shift)
                total += np.sum(((prediction[valid] - observed["value"][valid]) /
                                 observed["sigma"][valid]) ** 2)
            chi2[gi, si] = total
    best_index = np.unravel_index(np.argmin(chi2), chi2.shape)
    return model, shifts, chi2, float(gammas[best_index[0]]), float(shifts[best_index[1]])


def display_models(mesh, datasets, roughness_fraction, model_phases=64,
                   crater_theta_bins=8, crater_azimuth_bins=16,
                   facet_chunk_size=64, self_heating_iterations=4,
                   global_shadowing=True, global_self_heating=True,
                   view_factor_cache: Path | None = None):
    """独立计算论文图示的 Γ=330、350、370 三条展示曲线。"""
    output = {day: {} for day in datasets}
    base = ThermoConfig(n_phase=model_phases, crater_theta_bins=crater_theta_bins,
                        crater_azimuth_bins=crater_azimuth_bins,
                        facet_chunk_size=facet_chunk_size,
                        crater_self_heating_iterations=self_heating_iterations,
                        global_shadowing=global_shadowing,
                        global_self_heating=global_self_heating,
                        global_view_factor_cache=(str(view_factor_cache)
                                                  if view_factor_cache else None))
    for day in datasets:
        config = replace(base, **DATE_CONFIG[day])
        geometry = geometry_on_model_grid(datasets[day], model_phases)
        for gamma in (330.0, 350.0, 370.0):
            print(f"Figure model: {day}, Gamma={gamma:g}", flush=True)
            result = (simulate(mesh, gamma, config) if geometry is None else
                      simulate(mesh, gamma, config, *geometry))
            curve4 = result.mixed_detector(roughness_fraction)
            output[day][gamma] = curve4 / 1.0e4
    return output


def write_binned_csv(path, datasets):
    """写出可审计的一度分箱观测、误差、填充率及几何表。"""
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["date", "rotation_phase",
                         "detector_thermal_radiance_W_cm-2_um-1_sr-1",
                         "sigma", "frame_count", "PDS_FILL_FAC",
                         "disk_equivalent_radiance_diagnostic_W_cm-2_um-1_sr-1",
                         "sun_x_IAU_BENNU", "sun_y_IAU_BENNU",
                         "sun_z_IAU_BENNU", "observer_x_IAU_BENNU",
                         "observer_y_IAU_BENNU", "observer_z_IAU_BENNU",
                         "heliocentric_distance_au"])
        for day, data in datasets.items():
            for row in zip(data["phase"], data["value"], data["sigma"], data["count"],
                           data["pds_fill_factor"], data["disk_equivalent_diagnostic"],
                           data["sun_direction"][:, 0], data["sun_direction"][:, 1],
                           data["sun_direction"][:, 2], data["observer_direction"][:, 0],
                           data["observer_direction"][:, 1], data["observer_direction"][:, 2],
                           data["heliocentric_distance_au"]):
                writer.writerow([day, *row])


def write_figure(path, datasets, model, best_shift):
    """不依赖绘图库，直接生成补充图 2b 风格的 SVG。"""
    width, height = 920, 640
    left, right, top, bottom = 115, 875, 65, 545
    plot_phase = np.linspace(0, 1, 721)
    observed_values = np.concatenate([x["value"][np.isfinite(x["value"])] for x in datasets.values()])
    ymin = min(1.0e-4, float(np.min(observed_values)) * 0.97)
    ymax = max(1.3e-4, float(np.max(observed_values)) * 1.03)
    sx = lambda x: left + x * (right - left)
    sy = lambda y: bottom - (y - ymin) / (ymax - ymin) * (bottom - top)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#111}.axis{stroke:#111;stroke-width:1.5}.grid{stroke:#ddd}</style>']
    for tick in np.linspace(0, 1, 6):
        x = sx(tick)
        parts += [f'<line class="grid" x1="{x}" y1="{top}" x2="{x}" y2="{bottom}"/>',
                  f'<text x="{x}" y="{bottom+28}" text-anchor="middle" font-size="17">{tick:.1f}</text>']
    for tick in np.arange(np.floor(ymin/0.05e-4)*0.05e-4, ymax+0.01e-4, 0.05e-4):
        y = sy(tick)
        parts += [f'<line class="grid" x1="{left}" y1="{y}" x2="{right}" y2="{y}"/>',
                  f'<text x="{left-12}" y="{y+6}" text-anchor="end" font-size="17">{tick/1e-4:.2f}</text>']
    parts += [f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}"/>',
              f'<line class="axis" x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}"/>',
              f'<text x="{(left+right)/2}" y="{height-35}" text-anchor="middle" font-size="20">Rotation phase</text>',
              f'<text x="28" y="{(top+bottom)/2}" transform="rotate(-90 28 {(top+bottom)/2})" text-anchor="middle" font-size="19">Radiance (10^-4 W cm^-2 um^-1 sr^-1)</text>',
              f'<text x="{right-10}" y="{top+32}" text-anchor="end" font-size="28">(b)</text>']
    colors = {330.0: ("#22cc33", "8 7"), 350.0: ("#e11", ""), 370.0: ("#111", "")}
    legend_y = top + 28
    parts.append(f'<circle cx="{left+17}" cy="{legend_y}" r="4" fill="#0645e5"/><text x="{left+38}" y="{legend_y+6}" font-size="18">OVIRS 4 um data</text>')
    for row, gamma in enumerate((330.0, 350.0, 370.0), 1):
        color, dash = colors[gamma]
        curves = [cyclic_interpolate(model[day][gamma], plot_phase - best_shift) for day in datasets]
        curve = np.mean(curves, axis=0)
        points = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in zip(plot_phase, curve))
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.5"{dash_attr}/>')
        y = legend_y + row * 28
        parts.append(f'<line x1="{left+2}" y1="{y}" x2="{left+32}" y2="{y}" stroke="{color}" stroke-width="2.5"{dash_attr}/><text x="{left+38}" y="{y+6}" font-size="18">{int(gamma)}</text>')
    for data in datasets.values():
        valid = np.isfinite(data["value"])
        parts.extend(f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="2.0" fill="#0645e5"/>'
                     for x, y in zip(data["phase"][valid], data["value"][valid]))
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def main():
    """解析命令行、执行数据处理—正演网格—χ² 反演—结果输出全流程。"""
    # 避免 Windows 默认 GBK 在帮助或进度信息中遇到 µ、χ² 时编码失败。
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/ovirs/full"))
    parser.add_argument("--shape", type=Path,
                        default=Path("data/shape/g_12560mm_spc_obj_0000n00000_v014.obj"),
                        help="Bennu OBJ (default: official SPCv14 12.56 m)")
    parser.add_argument("--output", type=Path, default=Path("output/supplementary_fig2b"))
    parser.add_argument("--gamma-min", type=float, default=0.0)
    parser.add_argument("--gamma-max", type=float, default=600.0)
    parser.add_argument("--gamma-step", type=float, default=10.0)
    parser.add_argument("--model-phases", type=int, default=64)
    parser.add_argument("--fit-facet-stride", type=int, default=16,
                        help="uniform plate subsampling for the full Gamma scan")
    parser.add_argument("--facet-chunk-size", type=int, default=64,
                        help="macro facets per temporary radiative-exchange block")
    parser.add_argument("--self-heating-iterations", type=int, default=4)
    parser.add_argument("--crater-theta-bins", type=int, default=8)
    parser.add_argument("--crater-azimuth-bins", type=int, default=16)
    parser.add_argument("--spice-kernels", type=Path,
                        default=Path(r"D:\ATPM\data\spice\kernels"),
                        help="ASCII-only path required by CSPICE on Windows")
    parser.add_argument("--disable-global-shadowing", action="store_true")
    parser.add_argument("--disable-global-self-heating", action="store_true")
    parser.add_argument("--view-factor-cache", type=Path,
                        default=Path("data/shape/global_view_factors_spcv14.npz"))
    parser.add_argument("--skip-full-shape-figure", action="store_true",
                        help="plot the fast scan mesh instead of recomputing 330/350/370 on all plates")
    parser.add_argument("--skip-figure", action="store_true",
                        help="skip figure rendering (useful for a restricted Gamma scan)")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows = process_frames(args.input, args.output / "processed_ovirs_frames.csv")
    rows = attach_spice_geometry(rows, args.spice_kernels,
                                 args.output / "processed_ovirs_frames.csv")
    datasets, epoch = bin_one_degree(rows)
    write_binned_csv(args.output / "ovirs_thermal_1deg.csv", datasets)
    mesh = load_obj(args.shape) if args.shape else make_bennu_proxy_mesh()
    fit_mesh = sampled_mesh(mesh, args.fit_facet_stride)
    gammas = np.arange(args.gamma_min, args.gamma_max + args.gamma_step / 2, args.gamma_step)
    model, shifts, chi2, best_gamma, best_shift = fit_gamma(
        fit_mesh, datasets, gammas, roughness_fraction=0.77,
        model_phases=args.model_phases, facet_chunk_size=args.facet_chunk_size,
        self_heating_iterations=args.self_heating_iterations,
        crater_theta_bins=args.crater_theta_bins,
        crater_azimuth_bins=args.crater_azimuth_bins,
        global_shadowing=not args.disable_global_shadowing,
        global_self_heating=not args.disable_global_self_heating,
        view_factor_cache=args.view_factor_cache
    )
    with (args.output / "ovirs_gamma_phase_chi2.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["thermal_inertia", "phase_shift", "chi2"])
        for gi, gamma in enumerate(gammas):
            writer.writerows((gamma, shift, chi2[gi, si]) for si, shift in enumerate(shifts))
    profile = np.min(chi2, axis=1)
    inside = gammas[profile <= np.min(profile) + 1.0]
    data_count = int(sum(np.count_nonzero(np.isfinite(data["value"]) &
                                         np.isfinite(data["sigma"]))
                         for data in datasets.values()))
    degrees_of_freedom = max(1, data_count - 2)
    reduced_chi2 = float(np.min(chi2) / degrees_of_freedom)
    scaled_inside = gammas[profile <= np.min(profile) + reduced_chi2]
    summary = {
        "status": "public PDS L2 reconstructed to paper-like OVIRS thermal series",
        "strict_paper_reproduction": False,
        "input_fits": len(list(args.input.glob("*.fits"))),
        "accepted_frames": len(rows),
        "binning": "1 degree rotational phase; fit uses reflected-light-subtracted L2 detector radiance without division by PDS FILL_FAC",
        "reflection_removal": "PDS OSIRIS-REx project solar spectrum scaled over 2.05-2.15 um; flat reflectance; Figure 2b samples averaged over 3.98-4.02 um",
        "radiance_products": {
            "detector": "public L2 radiance after reflected-sunlight subtraction",
            "fit_quantity": "detector thermal radiance, compared with the CK/IK aperture-integrated model",
            "disk_equivalent_diagnostic": "detector thermal radiance divided by PDS FILL_FAC; retained for diagnostics and not fitted",
            "underfilled_fov_warning": "public calibration documentation reports possible artifacts; no author L3a empirical correction array is public"
        },
        "roughness": {"crater_fraction": 0.77, "rms_slope_deg": 43.0},
        "roughness_physics": {
            "crater_aperture_ray_visibility": True,
            "view_factor_thermal_self_heating": True,
            "multiple_scattered_sunlight": True,
            "directional_thermal_beaming": True,
            "crater_wall_elements": int(args.crater_theta_bins * args.crater_azimuth_bins),
            "self_heating_iterations": int(args.self_heating_iterations),
            "global_shape_shadowing": not args.disable_global_shadowing,
            "global_shape_mutual_heating_smooth_endpoint": not args.disable_global_self_heating,
            "rough_global_coupling": (
                "smooth-endpoint external field approximation"
                if not args.disable_global_self_heating else "disabled")
        },
        "spice_geometry": {
            "enabled": True,
            "evaluated_frames": len(rows),
            "body_fixed_frame": "IAU_BENNU",
            "aberration_correction": "LT+S",
            "kernel_set": KERNEL_SET_ID,
            "free_rotation_phase_shift": False
        },
        "shape": mesh.source,
        "shape_facets": int(len(mesh.faces)),
        "fit_grid_facets": int(len(fit_mesh.faces)),
        "fit_facet_stride": int(args.fit_facet_stride),
        "best_fit": {"thermal_inertia": best_gamma, "phase_shift": best_shift,
                     "chi2": float(np.min(chi2)),
                     "degrees_of_freedom": degrees_of_freedom,
                     "reduced_chi2": reduced_chi2,
                     "delta_chi2_1_interval": ([float(inside.min()), float(inside.max())]
                                               if len(inside) else None),
                     "model_discrepancy_scaled_interval": (
                         [float(scaled_inside.min()), float(scaled_inside.max())]
                         if len(scaled_inside) else None)},
        "paper_reference": {"thermal_inertia": 350.0, "sigma": 20.0,
                            "figure_curves": [330.0, 350.0, 370.0]},
        "limitations": ["public archive contains L2 v2, not authors' processed L3a array",
                        "PDS FILL_FAC is a diagnostic occupancy fraction and is not divided into the fitted radiance",
                        "public calibration documentation notes underfilled-FOV artifacts but no author-specific L3a correction array is public",
                        "Gamma scan uniformly subsamples shape plates for runtime; plotted curves can use all plates",
                        "per-frame SPICE geometry is binned/interpolated to the thermal solver phase grid",
                        "global mutual heating is solved on the smooth macro-facet endpoint and reused as an external field for rough microfacets, not fully coupled",
                        "reduced chi-square above unity means the formal Delta-chi-square interval is not a reliable physical uncertainty"],
        "phase_epoch_utc": datetime.fromtimestamp(epoch, tz=timezone.utc).isoformat(),
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    # 图中的三条论文参考曲线与本次 Γ 扫描范围分开计算，因此窄范围审计扫描也能绘图。
    if not args.skip_figure:
        figure_mesh = fit_mesh if args.skip_full_shape_figure else mesh
        figure_model = display_models(
            figure_mesh, datasets, roughness_fraction=0.77, model_phases=args.model_phases,
            facet_chunk_size=args.facet_chunk_size,
            self_heating_iterations=args.self_heating_iterations,
            crater_theta_bins=args.crater_theta_bins,
            crater_azimuth_bins=args.crater_azimuth_bins,
            global_shadowing=not args.disable_global_shadowing,
            global_self_heating=not args.disable_global_self_heating,
            view_factor_cache=args.view_factor_cache)
        write_figure(args.output / "supplementary_figure_2b_reproduction.svg",
                     datasets, figure_model, best_shift)
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
