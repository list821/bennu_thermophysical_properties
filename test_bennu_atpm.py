import unittest
import numpy as np

from bennu_atpm import (ModelCurves, ThermoConfig, crater_fraction_to_rms_deg,
                        make_bennu_proxy_mesh, planck_lambda, planck_ovirs_response,
                        rms_deg_to_crater_fraction, simulate,
                        _crater_aperture_visibility, _crater_local_geometry,
                        _facet_bases, _hemispherical_crater_geometry,
                        _multiple_scattered_solar)


class BennuAtpmTests(unittest.TestCase):
    """用物理守恒量约束 ATPM 关键逻辑，防止后续重构悄然改变结果。"""

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

    def test_hemisphere_acf_and_projected_area(self):
        """90° 半球坑的显式 ACF 必须还原单位宏观面元投影面积。"""
        cfg = ThermoConfig(crater_theta_bins=8, crater_azimuth_bins=16)
        crater = _hemispherical_crater_geometry(cfg)
        self.assertEqual(len(crater.normals), 128)
        self.assertAlmostEqual(crater.projected_area, np.pi, places=12)
        self.assertAlmostEqual(
            crater.area_conversion_factor*crater.projected_area, 1.0, places=12)
        # 半球实际面积/坑口面积为 2；z 投影面积/坑口面积为 1。
        self.assertAlmostEqual(float(np.sum(crater.area_weights)), 2.0, places=12)
        self.assertAlmostEqual(
            float(np.sum(crater.area_weights*crater.normals[:, 2])), 1.0, places=12)
        # 理想半球坑壁看到其它坑壁的总份额趋近 1/2；离散时扣除自身单元。
        expected_row_sum = .5-.5/len(crater.normals)
        self.assertTrue(np.allclose(crater.view_factors.sum(axis=1),
                                    expected_row_sum, atol=1e-14))

    def test_global_to_local_basis_is_orthonormal(self):
        """宏观面元到粗糙微面元的局部坐标映射必须保持正交和手性。"""
        normals = np.array([[0., 0., 1.], [1., 2., 3.], [-2., 1., .2]])
        normals /= np.linalg.norm(normals, axis=1)[:, None]
        t1, t2 = _facet_bases(normals)
        self.assertTrue(np.allclose(np.sum(t1*normals, axis=1), 0.0, atol=1e-14))
        self.assertTrue(np.allclose(np.sum(t2*normals, axis=1), 0.0, atol=1e-14))
        self.assertTrue(np.allclose(np.sum(t1*t2, axis=1), 0.0, atol=1e-14))
        self.assertTrue(np.allclose(np.cross(t1, t2), normals, atol=1e-14))

    def test_equation_24_smooth_rough_mixing(self):
        """式(24)应在两个端元间线性插值，且端点不能被额外 ACF 缩放。"""
        phase = np.array([0., .5])
        curves = ModelCurves(phase, np.array([1., 2.]), np.array([3., 6.]),
                             np.array([10., 20.]), np.array([30., 60.]),
                             np.ones(2), 0.0)
        self.assertTrue(np.array_equal(curves.mixed(0.)[0], curves.smooth_4um))
        self.assertTrue(np.array_equal(curves.mixed(1.)[0], curves.crater_4um))
        self.assertTrue(np.allclose(curves.mixed(.5)[0], np.array([2., 4.])))


if __name__ == "__main__":
    unittest.main()
