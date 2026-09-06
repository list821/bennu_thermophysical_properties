import unittest
import numpy as np

from bennu_atpm import (ThermoConfig, crater_fraction_to_rms_deg,
                        make_bennu_proxy_mesh, planck_lambda, planck_ovirs_response,
                        rms_deg_to_crater_fraction, simulate,
                        _crater_aperture_visibility, _crater_local_geometry,
                        _multiple_scattered_solar)


class BennuAtpmTests(unittest.TestCase):
    def test_planck_monotonic(self):
        values = planck_lambda(4., np.array([200., 300., 400.]))
        self.assertTrue(np.all(values > 0))
        self.assertTrue(np.all(np.diff(values) > 0))
        temperature = np.array([[250.0, 300.0, 350.0]])
        convolved = planck_ovirs_response(temperature, 4.00038, 0.00641)
        self.assertEqual(convolved.shape, temperature.shape)
        self.assertTrue(np.all(np.diff(convolved[0]) > 0))

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

    def test_crater_view_factors_scattering_and_ray_visibility(self):
        cfg = ThermoConfig(crater_theta_bins=3, crater_azimuth_bins=6)
        positions, _, weights, view = _crater_local_geometry(cfg)
        self.assertTrue(np.allclose(np.diag(view), 0.0))
        self.assertTrue(np.all(view.sum(axis=1) < 0.5))
        self.assertTrue(np.all(view.sum(axis=1) > 0.0))
        direct = np.ones((1, len(weights), 2))
        absorbed_single = (1-0.2)*direct
        absorbed_multiple = _multiple_scattered_solar(direct, view, 0.2)
        self.assertGreater(float(absorbed_multiple.sum()),
                           float(absorbed_single.sum()))
        zenith = np.array([[[0.0, 0.0, 1.0]]])
        below = np.array([[[0.0, 0.0, -1.0]]])
        self.assertTrue(np.all(_crater_aperture_visibility(positions, zenith)))
        self.assertFalse(np.any(_crater_aperture_visibility(positions, below)))


if __name__ == "__main__":
    unittest.main()
