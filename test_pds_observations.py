import os
import unittest
from pathlib import Path

import numpy as np

from pds_observations import read_ovirs_product, extract_ovirs_thermal_band


ROOT = Path(__file__).resolve().parent
OVIRS_FILENAME = "20181102T041236S859_ovr_scil2_calv2.fits"


def ovirs_test_product() -> Path:
    """Locate the OVIRS fixture using an optional external data root."""
    data_root = Path(os.environ.get("OVIRS_DATA_ROOT", ROOT / "data" / "ovirs"))
    candidates = (
        data_root / "sample" / OVIRS_FILENAME,
        data_root / "full" / OVIRS_FILENAME,
        data_root / OVIRS_FILENAME,
    )
    for path in candidates:
        if path.is_file():
            return path
    searched = "\n  ".join(str(path) for path in candidates)
    raise unittest.SkipTest(
        "OVIRS test product was not found. Set OVIRS_DATA_ROOT to the OVIRS "
        f"data directory. Searched:\n  {searched}"
    )


class PDSObservationTests(unittest.TestCase):
    def test_ovirs_fits_layout(self):
        path = ovirs_test_product()
        product = read_ovirs_product(path)
        self.assertEqual(product["radiance"].shape, (20, 512))
        self.assertAlmostEqual(float(product["phase_angle"]), 5.15778448548787)
        thermal, sigma, lines = extract_ovirs_thermal_band(product)
        self.assertGreater(lines, 0)
        self.assertTrue(np.isfinite(thermal))
        self.assertTrue(np.isfinite(sigma))


if __name__ == "__main__":
    unittest.main()
