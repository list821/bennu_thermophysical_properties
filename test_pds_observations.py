import unittest
from pathlib import Path

import numpy as np

from pds_observations import read_ovirs_product, extract_ovirs_thermal_band


ROOT = Path(__file__).resolve().parent


class PDSObservationTests(unittest.TestCase):
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
