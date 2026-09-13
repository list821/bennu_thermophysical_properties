import unittest
import numpy as np

from reproduce_supplementary_fig2b import residual_fit_metrics
from bennu_atpm import (ModelCurves, ThermoConfig, crater_fraction_to_rms_deg,
                        make_bennu_proxy_mesh, planck_lambda, planck_ovirs_response,
                        ovirs_aperture_fill, rms_deg_to_crater_fraction, simulate,
                        _crater_aperture_visibility, _crater_local_geometry,
                        _facet_bases, _hemispherical_crater_geometry,
                        _multiple_scattered_solar, _periodic_temperature, SIGMA,
                        SOLAR_CONSTANT)


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

    def test_mae_selection_does_not_use_variance(self):
        """默认 MAE 只看绝对残差；改变 sigma 不得改变 MAE/RMSE。"""
        prediction = np.array([1.0, 3.0])
        observation = np.array([2.0, 1.0])
        first = residual_fit_metrics(
            prediction, observation, np.array([0.1, 10.0]))
        second = residual_fit_metrics(
            prediction, observation, np.array([10.0, 0.1]))
        self.assertAlmostEqual(first["mae"], 1.5)
        self.assertAlmostEqual(first["rmse"], np.sqrt(2.5))
        self.assertEqual(first["mae"], second["mae"])
        self.assertEqual(first["rmse"], second["rmse"])
        self.assertNotEqual(first["chi2"], second["chi2"])

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
        self.assertTrue(curves.thermal_solver_converged)

    def test_periodic_solver_reports_final_residual(self):
        cfg = ThermoConfig(n_phase=96)
        phase = np.arange(cfg.n_phase)/cfg.n_phase
        absorbed = (900*np.maximum(np.cos(2*np.pi*phase), 0))[None, :]
        solution = _periodic_temperature(
            absorbed, 350.0, cfg, max_iter=1, tolerance=1e-12)
        omega = 2*np.pi*np.fft.rfftfreq(
            cfg.n_phase, d=cfg.period_s/cfg.n_phase)
        conduction = np.fft.irfft(
            np.fft.rfft(solution.temperature, axis=1) *
            (350.0*np.sqrt(1j*omega))[None, :],
            n=cfg.n_phase, axis=1)
        actual = np.max(np.abs(
            absorbed-cfg.emissivity*SIGMA*solution.temperature**4-conduction))
        self.assertEqual(solution.iterations, 1)
        self.assertFalse(solution.converged)
        self.assertAlmostEqual(solution.residual_max_w_m2, float(actual), places=12)

    def test_periodic_solver_harmonic_admittance_sign_and_period(self):
        cfg = ThermoConfig(n_phase=384)
        phase = np.arange(cfg.n_phase)/cfg.n_phase
        omega = 2*np.pi/cfg.period_s
        expected = 250.0+20.0*np.cos(2*np.pi*phase)
        conduction = (350.0*20.0*np.sqrt(omega/2) *
                      (np.cos(2*np.pi*phase)-np.sin(2*np.pi*phase)))
        absorbed = (cfg.emissivity*SIGMA*expected**4+conduction)[None, :]
        solution = _periodic_temperature(absorbed, 350.0, cfg)
        self.assertTrue(solution.converged)
        self.assertLess(np.max(np.abs(solution.temperature[0]-expected)), 0.01)

    def test_phase_resolution_changes_fit_not_boundary_residual(self):
        """Sharp shadows expose grid error even when every solve has a small residual."""
        phase_counts = (96, 384, 650, 768)
        gamma_grid = np.arange(330.0, 371.0, 10.0)

        def forcing(n_phase):
            cfg = ThermoConfig(n_phase=n_phase)
            phase = np.arange(n_phase)/n_phase
            sunlight = np.maximum(np.cos(2*np.pi*phase), 0)
            visible = ~(((phase > .91) | (phase < .035)) |
                        ((phase > .105) & (phase < .155)))
            return cfg, ((1-cfg.bond_albedo)*SOLAR_CONSTANT /
                         cfg.heliocentric_distance_au**2 *
                         sunlight*visible)[None, :]

        reference_cfg, reference_q = forcing(768)
        reference = _periodic_temperature(
            reference_q, 350.0, reference_cfg)
        self.assertTrue(reference.converged)
        reference_phase = np.arange(769)/768
        reference_radiance = planck_lambda(
            4.00038, np.r_[reference.temperature[0], reference.temperature[0, 0]])
        best = {}
        for n_phase in phase_counts:
            cfg, absorbed = forcing(n_phase)
            phase = np.arange(n_phase)/n_phase
            observation = np.interp(phase, reference_phase, reference_radiance)
            scores = []
            for gamma in gamma_grid:
                solution = _periodic_temperature(absorbed, gamma, cfg)
                self.assertTrue(solution.converged)
                prediction = planck_lambda(4.00038, solution.temperature[0])
                scores.append(np.mean((prediction/observation-1)**2))
            best[n_phase] = float(gamma_grid[np.argmin(scores)])
        self.assertNotEqual(best[96], 350.0)
        self.assertEqual(best[384], 350.0)
        self.assertEqual(best[650], 350.0)
        self.assertEqual(best[768], 350.0)

    def test_ovirs_aperture_fill_decreases_off_axis(self):
        """CK/IK 孔径耦合在中心指向最大，偏到响应边缘后必须降低。"""
        mesh = make_bennu_proxy_mesh(8, 16)
        observer = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        centered = np.array([[-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
        offset = 1.5e-3
        boresight = centered.copy()
        boresight[1] = [-np.cos(offset), np.sin(offset), 0.0]
        fill = ovirs_aperture_fill(
            mesh, observer, np.array([190.0, 190.0]), boresight,
            np.array([2.01e-3, 2.01e-3]))
        self.assertGreater(fill[0], 0.0)
        self.assertGreater(fill[0], fill[1])

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
