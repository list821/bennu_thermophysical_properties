import unittest
import numpy as np

from bennu_atpm import (ThermoConfig, crater_fraction_to_rms_deg,
                        make_bennu_proxy_mesh, planck_lambda,
                        rms_deg_to_crater_fraction, simulate)


class BennuAtpmTests(unittest.TestCase):
    def test_planck_monotonic(self):
        values = planck_lambda(4., np.array([200., 300., 400.]))
        self.assertTrue(np.all(values > 0))
        self.assertTrue(np.all(np.diff(values) > 0))

    def test_paper_roughness_mapping(self):
        self.assertAlmostEqual(crater_fraction_to_rms_deg(.77), 43., places=10)
        self.assertAlmostEqual(rms_deg_to_crater_fraction(43.), .77, places=10)

    def test_mesh_and_forward_model(self):
        mesh = make_bennu_proxy_mesh(5, 10)
        self.assertTrue(np.all(mesh.areas > 0))
        cfg = ThermoConfig(n_phase=32, crater_theta_bins=2, crater_azimuth_bins=4)
        curves = simulate(mesh, 350., cfg)
        self.assertEqual(curves.smooth_4um.shape, (32,))
        self.assertTrue(np.all(np.isfinite(curves.crater_14um)))
        self.assertLess(curves.max_boundary_residual_w_m2, 1.)


if __name__ == "__main__":
    unittest.main()
