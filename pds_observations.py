"""读取 OSIRIS-REx 公开 OTES/OVIRS PDS4 数据并构建光变。

二进制偏移和数据类型遵循配套 PDS4 XML 标签，所有观测量均来自归档产品，
不从论文图片反向取点。除 NumPy 外只依赖 Python 标准库。
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np


OTES_DTYPE = np.dtype(
    {
        "names": ["sclk", "sclk_sub", "ick", "quality", "radiance",
                  "bt_uncertainty", "max_brightness_temperature", "wavenumber"],
        "formats": ["<u4", "<u2", "<u2", "<u2", ("<f4", 349),
                    "<f4", "<f4", ("<f4", 349)],
        "offsets": [0, 4, 6, 8, 10, 1406, 1410, 1414],
        "itemsize": 2810,
    }
)


@dataclass
class BinnedLightcurve:
    """按自转相位分箱后的单仪器光变和数据质量统计。"""
    instrument: str
    source: str
    wavelength_um: float
    phase: np.ndarray
    value: np.ndarray
    sigma: np.ndarray
    count: np.ndarray
    raw_count: int
    accepted_count: int
    units_before_normalization: str


def _local(tag: str) -> str:
    """去掉 XML 命名空间，只保留本地标签名。"""
    return tag.rsplit("}", 1)[-1]


def _first_text(root: ET.Element, name: str, default: str | None = None) -> str | None:
    """返回 XML 中首个同名元素文本，缺失时返回默认值。"""
    for element in root.iter():
        if _local(element.tag) == name and element.text:
            return element.text.strip()
    return default


def _parse_utc(text: str) -> datetime:
    """将 PDS UTC 文本统一转换为带时区的 UTC datetime。"""
    parsed = datetime.fromisoformat(text.strip().replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _parse_sclk_string(text: str) -> tuple[int, int]:
    """拆分航天器时钟的整秒与亚秒计数。"""
    match = re.search(r"(?:\d+/)?(\d+)(?:\.(\d+))?", text.strip())
    if not match:
        raise ValueError(f"Unrecognized SCLK string: {text!r}")
    return int(match.group(1)), int(match.group(2) or 0)


def read_otes_calibrated(path: str | Path) -> np.ndarray:
    """按 OTES SIS 固定记录布局以内存映射读取标定光谱。"""
    path = Path(path)
    if path.stat().st_size % OTES_DTYPE.itemsize:
        raise ValueError(f"OTES file is not a whole number of 2810-byte records: {path}")
    # 只读内存映射避免把每段约 24 MB 序列完整装入内存，也适合页面文件较小的机器。
    return np.memmap(path, dtype=OTES_DTYPE, mode="r")


def _read_be_float(record: bytes, offset: int) -> float:
    """从记录指定偏移读取一个大端 32 位浮点数。"""
    return float(np.frombuffer(record[offset:offset + 4], dtype=">f4", count=1)[0])


def _read_be_double(record: bytes, offset: int) -> float:
    """从记录指定偏移读取一个大端 64 位浮点数。"""
    return float(np.frombuffer(record[offset:offset + 8], dtype=">f8", count=1)[0])


def read_otes_geometry(label_path: str | Path) -> list[dict[str, object]]:
    """依据 PDS4 XML 描述读取 OTES FITS 内嵌几何二进制表。"""
    label_path = Path(label_path)
    root = ET.parse(label_path).getroot()
    table = next(e for e in root.iter() if _local(e.tag) == "Table_Binary")
    offset = int(next(e.text for e in table if _local(e.tag) == "offset"))
    records = int(next(e.text for e in table if _local(e.tag) == "records"))
    record_binary = next(e for e in table if _local(e.tag) == "Record_Binary")
    length = int(next(e.text for e in record_binary if _local(e.tag) == "record_length"))
    fits_name = _first_text(root, "file_name")
    data = (label_path.parent / str(fits_name)).read_bytes()
    if offset + records * length > len(data):
        raise ValueError(f"Geometry table exceeds FITS file: {label_path}")
    output = []
    for index in range(records):
        row = data[offset + index * length: offset + (index + 1) * length]
        sclk_text = row[0:51].decode("ascii", errors="ignore").strip(" \x00")
        utc_text = row[51:102].decode("ascii", errors="ignore").strip(" \x00")
        look_type = row[227:247].decode("ascii", errors="ignore").strip(" \x00")
        sclk, sub = _parse_sclk_string(sclk_text)
        output.append(
            {
                "sclk": sclk,
                "sclk_sub": sub,
                "utc": utc_text,
                "latitude_deg": _read_be_float(row, 102),
                "longitude_deg": _read_be_float(row, 106),
                "bore_flag": int(np.frombuffer(row[193:195], dtype=">i2", count=1)[0]),
                "incidence_deg": _read_be_float(row, 195),
                "emission_deg": _read_be_float(row, 199),
                "phase_deg": _read_be_float(row, 203),
                "target_range_km": _read_be_double(row, 211),
                "fov_diameter_distance": _read_be_float(row, 223),
                "look_type": look_type,
            }
        )
    return output


def load_otes_geometry_directory(path: str | Path) -> dict[tuple[int, int], dict[str, object]]:
    """合并目录内所有 OTES 几何表，并以航天器时钟作为索引。"""
    result: dict[tuple[int, int], dict[str, object]] = {}
    for label in sorted(Path(path).glob("*_ote_geo.xml")):
        for record in read_otes_geometry(label):
            result[(int(record["sclk"]), int(record["sclk_sub"]))] = record
    return result


def _nearest_geometry(
    geometry: dict[tuple[int, int], dict[str, object]],
    by_second: dict[int, list[tuple[int, dict[str, object]]]],
    sclk: int,
    sub: int,
) -> dict[str, object] | None:
    """优先精确匹配 SCLK；必要时在同一秒内取最接近的亚秒记录。"""
    exact = geometry.get((int(sclk), int(sub)))
    if exact is not None:
        return exact
    # 部分标签对 SCLK 小数舍入不同；回退匹配严格限制在同一整秒。
    candidates = [(abs(candidate_sub - int(sub)), value)
                  for candidate_sub, value in by_second.get(int(sclk), [])]
    return min(candidates, default=(0, None), key=lambda item: item[0])[1]


def _bin_normalized_curve(
    phase: np.ndarray,
    value: np.ndarray,
    n_bins: int,
    systematic_floor: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """把光变归一化到均值 1，并输出分箱均值、稳健误差和计数。"""
    index = np.floor(np.mod(phase, 1.0) * n_bins).astype(int)
    centers = (np.arange(n_bins) + 0.5) / n_bins
    mean = np.full(n_bins, np.nan)
    sigma = np.full(n_bins, np.nan)
    count = np.zeros(n_bins, int)
    scale = float(np.nanmean(value))
    normalized = value / scale
    for i in range(n_bins):
        samples = normalized[index == i]
        samples = samples[np.isfinite(samples)]
        count[i] = len(samples)
        if len(samples):
            mean[i] = float(np.mean(samples))
            if len(samples) > 1:
                # 使用 MAD 稳健标准误，同时保留仪器/系统误差下限。
                mad = 1.4826 * np.median(np.abs(samples - np.median(samples)))
                sigma[i] = max(float(mad / np.sqrt(len(samples))), systematic_floor)
            else:
                sigma[i] = systematic_floor
    return centers, mean, sigma, count


def load_otes_lightcurves(
    science_dir: str | Path,
    geometry_dir: str | Path,
    period_hours: float = 4.296061,
    wavelength_um: float = 14.0,
    band_half_width_cm: float = 18.0,
    n_bins: int = 96,
) -> tuple[list[BinnedLightcurve], dict[str, object]]:
    """质量筛选 OTES 约 14 µm 数据并生成各观测序列的归一化光变。"""
    geometry = load_otes_geometry_directory(geometry_dir)
    geometry_by_second: dict[int, list[tuple[int, dict[str, object]]]] = {}
    for (second, sub), item in geometry.items():
        geometry_by_second.setdefault(second, []).append((sub, item))
    target_wavenumber = 10000.0 / wavelength_um
    curves: list[BinnedLightcurve] = []
    diagnostics: dict[str, object] = {"geometry_records": len(geometry), "files": []}
    period_s = period_hours * 3600.0
    for path in sorted(Path(science_dir).glob("*_ote_scil2.dat")):
        records = read_otes_calibrated(path)
        values, times, phases, matched = [], [], [], 0
        selected_geometry = []
        for record in records:
            geo = _nearest_geometry(
                geometry, geometry_by_second, int(record["sclk"]), int(record["sclk_sub"])
            )
            if geo is None:
                continue
            matched += 1
            quality = int(record["quality"])
            # 低两位是辐射质量（0 最好，3 表示无太空定标）；第 3 位标记相位反转/
            # 无效的亮温光谱。
            quality_ok = (quality & 0b11) <= 2 and (quality & 0b100) == 0
            if not quality_ok or geo["look_type"].lower() != "data-look" or int(geo["bore_flag"]) != 1:
                continue
            wn = np.asarray(record["wavenumber"], float)
            rad = np.asarray(record["radiance"], float)
            band = np.isfinite(wn) & np.isfinite(rad) & (rad > 0) & (np.abs(wn - target_wavenumber) <= band_half_width_cm)
            if np.count_nonzero(band) < 2:
                continue
            values.append(float(np.mean(rad[band])))
            times.append(float(record["sclk"]) + float(record["sclk_sub"]) / 65536.0)
            selected_geometry.append(geo)
        values_array = np.asarray(values)
        times_array = np.asarray(times)
        if len(values_array) < n_bins:
            raise ValueError(f"Too few valid OTES data looks in {path}: {len(values_array)}")
        phase_array = np.mod((times_array - times_array.min()) / period_s, 1.0)
        center, mean, sigma, count = _bin_normalized_curve(
            phase_array, values_array, n_bins=n_bins, systematic_floor=0.005
        )
        curves.append(
            BinnedLightcurve(
                instrument="OTES",
                source=path.name,
                wavelength_um=wavelength_um,
                phase=center,
                value=mean,
                sigma=sigma,
                count=count,
                raw_count=len(records),
                accepted_count=len(values_array),
                units_before_normalization="W cm-2 sr-1 (cm-1)-1",
            )
        )
        phases.extend(float(g["phase_deg"]) for g in selected_geometry)
        diagnostics["files"].append(
            {
                "file": path.name,
                "raw_records": len(records),
                "geometry_matches": matched,
                "accepted_data_looks": len(values_array),
                "filled_phase_bins": int(np.count_nonzero(count)),
                "mean_raw_band_radiance": float(np.mean(values_array)),
                "median_phase_angle_deg": float(np.median(phases)),
            }
        )
    return curves, diagnostics


def _read_array(path: Path, offset: int, dtype: str, shape: tuple[int, int]) -> np.ndarray:
    """从 FITS 固定字节偏移读取指定形状的原始数组。"""
    with path.open("rb") as stream:
        stream.seek(offset)
        result = np.fromfile(stream, dtype=np.dtype(dtype), count=shape[0] * shape[1])
    if result.size != shape[0] * shape[1]:
        raise ValueError(f"Short array read at byte {offset}: {path}")
    return result.reshape(shape)


def _fits_primary_header(path: Path) -> dict[str, object]:
    """解析 OVIRS 5760 字节主头中的简单标量 FITS 卡片。"""
    with path.open("rb") as stream:
        raw = stream.read(5760)
    header: dict[str, object] = {}
    for start in range(0, len(raw), 80):
        card = raw[start:start + 80].decode("ascii", errors="replace")
        key = card[:8].strip()
        if key == "END":
            break
        if not key or card[8:10] != "= ":
            continue
        token = card[10:].split("/", 1)[0].strip()
        if token.startswith("'"):
            value: object = token.strip("'").strip()
        elif token in ("T", "F"):
            value = token == "T"
        else:
            try:
                value = float(token) if any(x in token.upper() for x in (".", "E")) else int(token)
            except ValueError:
                value = token
        header[key] = value
    return header


def read_ovirs_product(path: str | Path) -> dict[str, object]:
    """读取 OVIRS L2 v2 的辐亮度、质量、波长、不确定度和观测元数据。"""
    path = Path(path)
    label = path.with_suffix(".xml")
    header = _fits_primary_header(path)
    # v2 Approach 产品使用文档规定的五个 HDU 固定布局。
    radiance = _read_array(path, 5760, ">f8", (20, 512))
    quality = _read_array(path, 92160, ">i4", (20, 512))
    wavelength = _read_array(path, 138240, ">f8", (20, 512))
    uncertainty = _read_array(path, 475200, ">f8", (20, 512))
    if label.exists():
        root = ET.parse(label).getroot()
        mid_obs = _first_text(root, "mid_obs") or _first_text(root, "start_date_time")
        spatial_names = ("bore_flag", "fov_fill_factor", "phase_angle", "sun_range", "target_range")
        spatial = {name: float(_first_text(root, name, "nan")) for name in spatial_names}
    else:
        mid_obs = str(header["MIDOBS"])
        spatial = {
            "bore_flag": float(header.get("BS_FLAG", -1)),
            "fov_fill_factor": float(header.get("FILL_FAC", np.nan)),
            "phase_angle": float(header.get("PHASEANG", np.nan)),
            "sun_range": float(header.get("SUN_RNG", np.nan)),
            "target_range": float(header.get("TARGRNG", np.nan)),
        }
    return {
        "path": str(path),
        "mid_obs": _parse_utc(str(mid_obs)),
        "radiance": radiance,
        "quality": quality,
        "wavelength_um": wavelength,
        "uncertainty": uncertainty,
        **spatial,
    }


def _planck_shape_um(wavelength_um: np.ndarray, temperature: float = 5778.0) -> np.ndarray:
    """没有归档太阳谱时使用的 5778 K Planck 相对光谱后备形状。"""
    c2_um_k = 1.438776877e4
    wavelength_um = np.asarray(wavelength_um, float)
    return 1.0 / (wavelength_um**5 * np.expm1(c2_um_k / (wavelength_um * temperature)))


def _project_solar_shape(wavelength_um: np.ndarray) -> np.ndarray:
    """把归档 OSIRIS-REx 项目太阳通量谱插值到各 OVIRS 波长。"""
    path = Path(__file__).resolve().parent / "data" / "ovirs" / "orexsolarflux.csv"
    if not path.exists():
        return _planck_shape_um(wavelength_um)
    table = np.loadtxt(path, delimiter=",")
    return np.interp(wavelength_um, table[:, 0], table[:, 1], left=np.nan, right=np.nan)


def extract_ovirs_thermal_band(product: dict[str, object]) -> tuple[float, float, int]:
    """返回扣除反射太阳光后的 OVIRS 4 µm 辐亮度。

    论文反演使用 3.5–4.0 µm 全谱，但补充图 2b 明确是 4 µm 光变。因此单波长
    正演应和最接近 4 µm 的样本比较，不能在 Wien 侧陡峭的整个区间直接平均。
    """
    radiance = np.asarray(product["radiance"], float)
    quality = np.asarray(product["quality"])
    wavelength = np.asarray(product["wavelength_um"], float)
    uncertainty = np.asarray(product["uncertainty"], float)
    # 低半字节给超像素中的良好探测器像元数；第 4 位为空超像素，第 6 位为被拒异常值。
    valid_quality = ((quality & 0x0F) > 0) & ((quality & 0x10) == 0) & ((quality & 0x40) == 0)
    valid = valid_quality & np.isfinite(radiance) & np.isfinite(wavelength)
    # 按论文在约 2.1 µm 缩放太阳谱；对称窄区间可避免比例被蓝侧偏置。
    anchor = valid & (wavelength >= 2.05) & (wavelength <= 2.15)
    thermal = valid & (wavelength >= 3.98) & (wavelength <= 4.02)
    if np.count_nonzero(anchor) < 4 or np.count_nonzero(thermal) < 4:
        return np.nan, np.nan, 0
    solar = _project_solar_shape(wavelength)
    valid_solar = np.isfinite(solar) & (solar > 0)
    anchor &= valid_solar
    thermal &= valid_solar
    weight = np.ones_like(uncertainty, dtype=float)
    usable_uncertainty = np.isfinite(uncertainty) & (uncertainty > 0)
    np.divide(1.0, uncertainty**2, out=weight, where=usable_uncertainty)
    scale = np.sum(weight[anchor] * radiance[anchor] * solar[anchor]) / np.sum(
        weight[anchor] * solar[anchor] ** 2
    )
    residual = radiance[thermal] - scale * solar[thermal]
    finite_unc = uncertainty[thermal][np.isfinite(uncertainty[thermal])]
    propagated = float(np.sqrt(np.sum(finite_unc**2)) / max(len(finite_unc), 1))
    scatter = float(np.nanstd(residual, ddof=1) / np.sqrt(len(residual))) if len(residual) > 1 else 0.0
    return float(np.mean(residual)), max(propagated, scatter), int(np.count_nonzero(thermal))


def load_ovirs_samples(path: str | Path) -> list[dict[str, object]]:
    """批量提取目录内每个 OVIRS FITS 的 4 µm 热辐亮度与诊断量。"""
    output = []
    for fits in sorted(Path(path).glob("*_ovr_scil2_calv2.fits")):
        product = read_ovirs_product(fits)
        thermal, sigma, lines = extract_ovirs_thermal_band(product)
        output.append(
            {
                "file": fits.name,
                "mid_obs": product["mid_obs"].isoformat(),
                "thermal_4um_radiance": thermal,
                "sigma": sigma,
                "accepted_superpixels": lines,
                "phase_angle_deg": product["phase_angle"],
                "fov_fill_factor": product["fov_fill_factor"],
                "sun_range_km": product["sun_range"],
                "target_range_km": product["target_range"],
            }
        )
    return output


def load_ovirs_lightcurves(
    path: str | Path,
    period_hours: float = 4.296061,
    n_bins: int = 96,
    minimum_samples: int | None = None,
) -> tuple[list[BinnedLightcurve], dict[str, object]]:
    """构建逐日的扣反射光 OVIRS 4 µm 光变。

    反射分量以 2.05–2.15 µm 太阳光谱形状拟合。PDS FILL_FAC 只作为观测元数据，
    不再充当额外辐射除数；保留标定的检测器辐亮度。样本不足的日期只报告不拟合。
    """
    minimum_samples = n_bins if minimum_samples is None else minimum_samples
    by_day: dict[str, list[tuple[datetime, float, float, str]]] = {}
    diagnostics: dict[str, object] = {"files_seen": 0, "files_accepted": 0, "days": []}
    for fits in sorted(Path(path).glob("*_ovr_scil2_calv2.fits")):
        diagnostics["files_seen"] += 1
        product = read_ovirs_product(fits)
        thermal, sigma, samples = extract_ovirs_thermal_band(product)
        fill = float(product["fov_fill_factor"])
        if (not np.isfinite(thermal) or thermal <= 0 or not np.isfinite(fill) or fill <= 0
                or samples == 0):
            continue
        # 补充图 2b 使用标定 OVIRS 辐亮度，论文未要求再除几何 FILL_FAC。
        corrected = thermal
        corrected_sigma = max(sigma, abs(corrected) * 0.01)
        timestamp = product["mid_obs"]
        day = timestamp.date().isoformat()
        by_day.setdefault(day, []).append((timestamp, corrected, corrected_sigma, fits.name))
        diagnostics["files_accepted"] += 1
    curves = []
    period_s = period_hours * 3600.0
    for day, rows in sorted(by_day.items()):
        rows.sort(key=lambda row: row[0])
        times = np.asarray([(row[0] - rows[0][0]).total_seconds() for row in rows])
        values = np.asarray([row[1] for row in rows])
        status = "used" if len(rows) >= minimum_samples else "insufficient samples"
        diagnostics["days"].append({"date": day, "accepted_samples": len(rows), "status": status})
        if len(rows) < minimum_samples:
            continue
        center, mean, empirical_sigma, count = _bin_normalized_curve(
            np.mod(times / period_s, 1.0), values, n_bins=n_bins, systematic_floor=0.01
        )
        curves.append(BinnedLightcurve(
            instrument="OVIRS", source=f"OVIRS {day}", wavelength_um=3.75,
            phase=center, value=mean, sigma=empirical_sigma, count=count,
            raw_count=len(rows), accepted_count=len(rows),
            units_before_normalization="W cm-2 sr-1 um-1 after FOV-fill correction",
        ))
    return curves, diagnostics
