"""使用 NAIF 官方 SPICE 核计算逐帧 OSIRIS-REx 观测几何。

方向向量分别从 Bennu 指向太阳和航天器，并表达在与 SPC 形状轴一致的
``IAU_BENNU`` 天体固定坐标系中。SPK 给位置，PCK/框架核给 Bennu 姿态，
CK 给航天器姿态，IK 给 OVIRS 实际视轴和圆形视场边缘。
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np

from download_spice_kernels import KERNELS

AU_KM = 149_597_870.7
BENNU_ID = "2101955"
SUN_ID = "10"
ORX_ID = "-64"
KERNEL_SET_ID = ("SPCv14-compatible: bennu_v14 PCK + reconstructed 2018 ORX SPK/CK "
                 "+ OVIRS v00 IK")


def load_kernel_set(kernel_root: str | Path):
    """检查并装载本项目固定的 SPICE 核集合，返回 spiceypy 接口。"""
    try:
        import spiceypy as spice
    except ImportError as exc:
        raise RuntimeError("SPICE geometry requires spiceypy") from exc
    root = Path(kernel_root)
    missing = [str(root / item) for item in KERNELS if not (root / item).exists()]
    if missing:
        raise FileNotFoundError("Missing SPICE kernels; run download_spice_kernels.py:\n" +
                                "\n".join(missing))
    spice.kclear()
    # 不调用 Path.resolve()：Windows 会把 ASCII 联接还原成含中文的项目路径，
    # 而 CSPICE N0067 对该路径的打开并不可靠。
    for relative in KERNELS:
        spice.furnsh(str(root.absolute() / relative))
    return spice


def frame_geometry(utc_values: Iterable[str], kernel_root: str | Path,
                   frame: str = "IAU_BENNU") -> dict[str, np.ndarray]:
    """在每个输入 UTC 时刻求太阳、观测者、OVIRS 视轴和视场半角。"""
    spice = load_kernel_set(kernel_root)
    utc = list(utc_values)
    # FITS 时间缓存为带 +00:00 的 ISO-8601；CSPICE 接收等价的无偏移或 Z 结尾 UTC。
    normalized_utc = [value[:-6] if value.endswith("+00:00") else value for value in utc]
    et = np.asarray([spice.str2et(value) for value in normalized_utc])
    try:
        # LT+S 同时修正单程光行时和恒星像差，与论文逐帧观测几何含义一致。
        sun_km, _ = spice.spkpos(SUN_ID, et, frame, "LT+S", BENNU_ID)
        observer_km, _ = spice.spkpos(ORX_ID, et, frame, "LT+S", BENNU_ID)
        shape, fov_frame, instrument_boresight, boundary_count, bounds = spice.getfov(-64321, 4)
        if shape != "CIRCLE" or boundary_count != 1:
            raise RuntimeError(f"Unexpected OVIRS IK FOV: {shape}, {boundary_count}")
        # 由 CK/框架链把 IK 中的仪器视轴逐帧旋转到 IAU_BENNU。
        transforms = np.asarray([spice.pxform(fov_frame, frame, value) for value in et])
        boresight = np.einsum("tij,j->ti", transforms, instrument_boresight)
        fov_half_angle = float(np.arccos(np.clip(
            np.dot(instrument_boresight, bounds[0]) /
            (np.linalg.norm(instrument_boresight)*np.linalg.norm(bounds[0])), -1, 1)))
    finally:
        spice.kclear()
    sun_km = np.asarray(sun_km, float)
    observer_km = np.asarray(observer_km, float)
    sun_distance = np.linalg.norm(sun_km, axis=1)
    observer_distance = np.linalg.norm(observer_km, axis=1)
    sun_direction = sun_km/sun_distance[:, None]
    observer_direction = observer_km/observer_distance[:, None]
    cosine = np.clip(np.einsum("ij,ij->i", sun_direction, observer_direction), -1, 1)
    return {
        "et_s": et,
        "sun_direction": sun_direction,
        "observer_direction": observer_direction,
        "spacecraft_position_km": observer_km,
        "ovirs_boresight_direction": boresight/np.linalg.norm(boresight, axis=1)[:, None],
        "ovirs_fov_half_angle_rad": np.full(len(et), fov_half_angle),
        "heliocentric_distance_au": sun_distance/AU_KM,
        "spacecraft_distance_km": observer_distance,
        "phase_angle_deg": np.rad2deg(np.arccos(cosine)),
        "subsolar_latitude_deg": np.rad2deg(np.arcsin(sun_direction[:, 2])),
    }
