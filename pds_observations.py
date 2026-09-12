"""读取 OSIRIS-REx 公开 OVIRS PDS4 数据并提取 4 µm 热辐亮度。

二进制偏移和数据类型遵循配套 PDS4 XML 标签，所有观测量均来自归档产品，
不从论文图片反向取点。除 NumPy 外只依赖 Python 标准库。
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np


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
            "fov_flag": float(header.get("FOV_FLAG", -1)),
            "boresight_angle_deg": float(header.get("BS_ANGLE", np.nan)),
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


def read_ovirs_pointing_metadata(path: str | Path) -> dict[str, float]:
    """只读 FITS 主头中的指向质量字段，不加载四个大型光谱数组。

    这个轻量入口用于给旧的逐帧光谱缓存补充 BS_FLAG、FOV_FLAG、
    BS_ANGLE 和 FILL_FAC，避免因增加质量筛选而重新处理全部光谱。
    """
    header = _fits_primary_header(Path(path))
    return {
        "bore_flag": float(header.get("BS_FLAG", -1)),
        "fov_flag": float(header.get("FOV_FLAG", -1)),
        "boresight_angle_deg": float(header.get("BS_ANGLE", np.nan)),
        "fov_fill_factor": float(header.get("FILL_FAC", np.nan)),
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
