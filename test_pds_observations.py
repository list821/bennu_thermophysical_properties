import unittest
from pathlib import Path

import numpy as np

from pds_observations import (load_otes_geometry_directory, read_otes_calibrated,
                              read_ovirs_product, extract_ovirs_thermal_band)


ROOT = Path(__file__).resolve().parent


class PDSObservationTests(unittest.TestCase):
    def test_otes_binary_layout(self):
        path = ROOT / "data/otes/science/20181108T040144S029_ote_scil2.dat"
        records = read_otes_calibrated(path)
        self.assertEqual(len(records), 8675)
        self.assertEqual(records["radiance"].shape, (8675, 349))
        self.assertTrue(np.isfinite(records["wavenumber"]).any())

    def test_otes_geometry(self):
        geometry = load_otes_geometry_directory(ROOT / "data/otes/geometry")
        self.assertGreater(len(geometry), 1000)
        self.assertTrue(any(str(item["look_type"]).lower() == "data-look" for item in geometry.values()))

    def test_ovirs_fits_layout(self):
        path = ROOT / "data/ovirs/sample/20181102T041236S859_ovr_scil2_calv2.fits"
        product = read_ovirs_product(path)
        self.assertEqual(product["radiance"].shape, (20, 512))
        self.assertAlmostEqual(float(product["phase_angle"]), 5.15778448548787)
        thermal, sigma, lines = extract_ovirs_thermal_band(product)
        self.assertGreater(lines, 0)
        self.assertTrue(np.isfinite(thermal))
        self.assertTrue(np.isfinite(sigma))


if __name__ == "__main__":
    unittest.main()
