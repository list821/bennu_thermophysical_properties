"""Readers and light-curve builders for public OSIRIS-REx OTES/OVIRS PDS4 data.

Only NumPy and the Python standard library are required.  Binary locations and
types follow the companion PDS4 XML labels; no values are scraped from plots.
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
    return tag.rsplit("}", 1)[-1]


def _first_text(root: ET.Element, name: str, default: str | None = None) -> str | None:
    for element in root.iter():
        if _local(element.tag) == name and element.text:
            return element.text.strip()
    return default


def _parse_utc(text: str) -> datetime:
    parsed = datetime.fromisoformat(text.strip().replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _parse_sclk_string(text: str) -> tuple[int, int]:
    match = re.search(r"(?:\d+/)?(\d+)(?:\.(\d+))?", text.strip())
    if not match:
        raise ValueError(f"Unrecognized SCLK string: {text!r}")
    return int(match.group(1)), int(match.group(2) or 0)


def read_otes_calibrated(path: str | Path) -> np.ndarray:
    path = Path(path)
    if path.stat().st_size % OTES_DTYPE.itemsize:
        raise ValueError(f"OTES file is not a whole number of 2810-byte records: {path}")
    # A read-only memory map avoids loading each ~24 MB observing sequence into
    # RAM and scales to machines with a small Windows page file.
    return np.memmap(path, dtype=OTES_DTYPE, mode="r")


def _read_be_float(record: bytes, offset: int) -> float:
    return float(np.frombuffer(record[offset:offset + 4], dtype=">f4", count=1)[0])


def _read_be_double(record: bytes, offset: int) -> float:
    return float(np.frombuffer(record[offset:offset + 8], dtype=">f8", count=1)[0])


def read_otes_geometry(label_path: str | Path) -> list[dict[str, object]]:
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
    exact = geometry.get((int(sclk), int(sub)))
    if exact is not None:
        return exact
    # Some labels round the fractional SCLK differently.  Restrict fallback to
    # the same second and choose the closest tick.
    candidates = [(abs(candidate_sub - int(sub)), value)
                  for candidate_sub, value in by_second.get(int(sclk), [])]
    return min(candidates, default=(0, None), key=lambda item: item[0])[1]


def _bin_normalized_curve(
    phase: np.ndarray,
    value: np.ndarray,
    n_bins: int,
    systematic_floor: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
                # Robust standard error; retain an instrumental/systematic floor.
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
            # Bits 1-2 are radiometric quality (0 best, 3 no space look); bit 3
            # marks a phase-inverted/invalid brightness-temperature spectrum.
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
    with path.open("rb") as stream:
        stream.seek(offset)
        result = np.fromfile(stream, dtype=np.dtype(dtype), count=shape[0] * shape[1])
    if result.size != shape[0] * shape[1]:
        raise ValueError(f"Short array read at byte {offset}: {path}")
    return result.reshape(shape)


def _fits_primary_header(path: Path) -> dict[str, object]:
    """Parse simple scalar cards from the 5760-byte OVIRS primary header."""
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
    path = Path(path)
    label = path.with_suffix(".xml")
    header = _fits_primary_header(path)
    # Version-2 approach products use the documented five-HDU fixed layout.
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
    c2_um_k = 1.438776877e4
    wavelength_um = np.asarray(wavelength_um, float)
    return 1.0 / (wavelength_um**5 * np.expm1(c2_um_k / (wavelength_um * temperature)))


def _project_solar_shape(wavelength_um: np.ndarray) -> np.ndarray:
    """Interpolate the archived OSIRIS-REx project solar-flux spectrum."""
    path = Path(__file__).resolve().parent / "data" / "ovirs" / "orexsolarflux.csv"
    if not path.exists():
        return _planck_shape_um(wavelength_um)
    table = np.loadtxt(path, delimiter=",")
    return np.interp(wavelength_um, table[:, 0], table[:, 1], left=np.nan, right=np.nan)


def extract_ovirs_thermal_band(product: dict[str, object]) -> tuple[float, float, int]:
    """Return the reflected-light-subtracted OVIRS radiance at 4 um.

    The paper fits spectra throughout 3.5--4.0 um, but Supplementary Fig. 2b
    is explicitly the 4-um light curve.  A single-wavelength model must
    therefore be compared with the samples nearest 4 um, not with an average
    across the steep Wien-side spectral interval.
    """
    radiance = np.asarray(product["radiance"], float)
    quality = np.asarray(product["quality"])
    wavelength = np.asarray(product["wavelength_um"], float)
    uncertainty = np.asarray(product["uncertainty"], float)
    # Low nibble is the count of good detector pixels in a superpixel.  Bit 4
    # marks an empty superpixel and bit 6 a rejected outlier (OVIRS SIS 4-4).
    valid_quality = ((quality & 0x0F) > 0) & ((quality & 0x10) == 0) & ((quality & 0x40) == 0)
    valid = valid_quality & np.isfinite(radiance) & np.isfinite(wavelength)
    # The paper scales the solar spectrum at approximately 2.1 um.  Use a
    # symmetric narrow interval to avoid biasing the scale to the blue side.
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
    """Build daily reflected-light-subtracted 4-um curves.

    The reflected component is fitted at 1.9-2.1 um with a solar spectral
    shape.  PDS FILL_FAC is treated as observing metadata, not as an extra
    radiometric divisor; the calibrated detector radiance is retained.
    One-file format samples are reported but intentionally not fitted.
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
        # Supplementary Fig. 2b uses the calibrated OVIRS radiance.  The paper
        # does not prescribe another division by the geometric FILL_FAC.
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
