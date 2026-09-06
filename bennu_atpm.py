"""Auditable, dependency-light thermophysical model for asteroid Bennu.

This implements triangular-facet illumination, periodic 1-D conduction,
fractional hemispherical-crater roughness, crater aperture ray visibility,
view-factor self-heating, multiple-scattered sunlight, Planck radiation and
directional disk integration.  It is an independent implementation of the
published Rozitis & Green (2011) method, not their unpublished C++ source.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

SIGMA = 5.670374419e-8
SOLAR_CONSTANT = 1361.0
H = 6.62607015e-34
C = 299792458.0
KB = 1.380649e-23


@dataclass(frozen=True)
class ThermoConfig:
    period_hours: float = 4.296061
    bond_albedo: float = 0.016
    emissivity: float = 0.9
    heliocentric_distance_au: float = 1.20
    density: float = 2500.0
    heat_capacity: float = 800.0
    phase_angle_deg: float = 5.0
    subsolar_latitude_deg: float = 0.0
    n_phase: int = 96
    # 128 equal-area wall elements is the production setting.  The former
    # 2x4 setting was useful for debugging but is too coarse for shadow-edge
    # and directional-beaming convergence.
    crater_theta_bins: int = 8
    crater_azimuth_bins: int = 16
    facet_chunk_size: int = 64
    crater_self_heating_iterations: int = 4
    crater_self_heating_tolerance_k: float = 0.08
    thermal_albedo: float = 0.0
    ovirs_wavelength_um: float = 4.00038
    ovirs_fwhm_um: float = 0.00641
    global_shadowing: bool = False
    global_self_heating: bool = False
    global_self_heating_iterations: int = 3
    global_view_factor_cache: str | None = None

    @property
    def period_s(self) -> float:
        return self.period_hours * 3600.0


@dataclass(frozen=True)
class Mesh:
    vertices: np.ndarray
    faces: np.ndarray
    centers: np.ndarray
    normals: np.ndarray
    areas: np.ndarray
    source: str


@dataclass
class ModelCurves:
    phase: np.ndarray
    smooth_4um: np.ndarray
    crater_4um: np.ndarray
    smooth_14um: np.ndarray
    crater_14um: np.ndarray
    projected_area_m2: np.ndarray
    max_boundary_residual_w_m2: float
    smooth_4um_detector: np.ndarray | None = None
    crater_4um_detector: np.ndarray | None = None
    modeled_fov_fill_factor: np.ndarray | None = None

    def mixed(self, crater_fraction: float) -> tuple[np.ndarray, np.ndarray]:
        f = float(np.clip(crater_fraction, 0.0, 1.0))
        return ((1-f)*self.smooth_4um + f*self.crater_4um,
                (1-f)*self.smooth_14um + f*self.crater_14um)

    def mixed_detector(self, crater_fraction: float) -> np.ndarray:
        if self.smooth_4um_detector is None or self.crater_4um_detector is None:
            return self.mixed(crater_fraction)[0]
        f = float(np.clip(crater_fraction, 0.0, 1.0))
        return (1-f)*self.smooth_4um_detector+f*self.crater_4um_detector


def _mesh_from_arrays(vertices: np.ndarray, faces: np.ndarray, source: str) -> Mesh:
    tri = vertices[faces]
    cross = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    centers = tri.mean(axis=1)
    inward = np.einsum("ij,ij->i", cross, centers) < 0
    if np.any(inward):
        faces = faces.copy()
        faces[inward] = faces[inward][:, [0, 2, 1]]
        tri = vertices[faces]
        cross = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
        centers = tri.mean(axis=1)
    norm = np.linalg.norm(cross, axis=1)
    good = norm > 1e-12
    return Mesh(vertices, faces[good], centers[good], cross[good]/norm[good, None],
                norm[good]/2, source)


def load_obj(path: str | Path) -> Mesh:
    vertices, faces = [], []
    obj_path = Path(path)
    for line in obj_path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.strip().split()
        if not fields:
            continue
        if fields[0] == "v" and len(fields) >= 4:
            vertices.append([float(x) for x in fields[1:4]])
        elif fields[0] == "f" and len(fields) >= 4:
            idx = [int(item.split("/")[0])-1 for item in fields[1:]]
            faces.extend([[idx[0], idx[j], idx[j+1]] for j in range(1, len(idx)-1)])
    if not vertices or not faces:
        raise ValueError(f"No usable vertices/faces in OBJ: {obj_path}")
    vertex_array = np.asarray(vertices, float)
    # PDS Bennu SPC OBJ products are explicitly in km; locally converted DSK
    # files and the proxy are in metres.  This conservative scale test keeps
    # projected areas physically meaningful without changing directions.
    scale_note = ""
    if np.max(np.linalg.norm(vertex_array, axis=1)) < 10.0:
        vertex_array *= 1000.0
        scale_note = " (PDS km coordinates converted to m)"
    return _mesh_from_arrays(vertex_array, np.asarray(faces, int),
                             str(obj_path.resolve()) + scale_note)


def make_bennu_proxy_mesh(n_latitude: int = 8, n_longitude: int = 16) -> Mesh:
    """Deterministic spinning-top proxy (~246 m radius), explicitly not SPC v13."""
    if n_latitude < 4 or n_longitude < 8:
        raise ValueError("Need at least 4 latitude and 8 longitude bins")
    vertices = [[0.0, 0.0, -220.0]]
    lats = np.linspace(-np.pi/2, np.pi/2, n_latitude+1)[1:-1]
    for lat in lats:
        for j in range(n_longitude):
            lon = 2*np.pi*j/n_longitude
            ridge = .105*np.cos(lat)**8
            flatten = -.055*np.sin(lat)**2
            irregular = (.018*np.sin(3*lon+.7)+.011*np.cos(5*lon-.4))*np.cos(lat)**2
            r = 246*(1+ridge+flatten+irregular)
            vertices.append([r*np.cos(lat)*np.cos(lon), r*np.cos(lat)*np.sin(lon), r*np.sin(lat)])
    north = len(vertices)
    vertices.append([0.0, 0.0, 220.0])
    faces, rings, first = [], len(lats), 1
    for j in range(n_longitude):
        faces.append([0, first+(j+1)%n_longitude, first+j])
    for ring in range(rings-1):
        a, b = 1+ring*n_longitude, 1+(ring+1)*n_longitude
        for j in range(n_longitude):
            k = (j+1)%n_longitude
            faces.extend([[a+j, a+k, b+k], [a+j, b+k, b+j]])
    last = 1+(rings-1)*n_longitude
    for j in range(n_longitude):
        faces.append([last+j, last+(j+1)%n_longitude, north])
    return _mesh_from_arrays(np.asarray(vertices), np.asarray(faces, int),
                             "internal Bennu spinning-top proxy (not SPC v13)")


def planck_lambda(wavelength_um: float | np.ndarray,
                  temperature_k: np.ndarray) -> np.ndarray:
    lam = np.asarray(wavelength_um, float)*1e-6
    temperature = np.maximum(np.asarray(temperature_k, float), 1.0)
    x = H*C/(lam*KB*temperature)
    return (2*H*C**2/lam**5/np.expm1(np.minimum(x, 700)))*1e-6


def planck_ovirs_response(temperature_k: np.ndarray, center_um: float,
                          fwhm_um: float, quadrature_order: int = 9) -> np.ndarray:
    """Gaussian OVIRS spectral-lineshape convolution at one spectral sample."""
    if fwhm_um <= 0:
        return planck_lambda(center_um, temperature_k)
    nodes, weights = np.polynomial.hermite.hermgauss(quadrature_order)
    sigma = fwhm_um/(2*np.sqrt(2*np.log(2)))
    wavelength = center_um+np.sqrt(2)*sigma*nodes
    temperature = np.asarray(temperature_k)
    wavelength_shape = (len(wavelength),)+(1,)*temperature.ndim
    values = planck_lambda(wavelength.reshape(wavelength_shape), temperature[None, ...])
    return np.einsum("q,q...->...", weights/np.sqrt(np.pi), values)


def skin_depth_m(gamma: float, config: ThermoConfig) -> float:
    return gamma/(config.density*config.heat_capacity)*np.sqrt(2/(2*np.pi/config.period_s))


def _periodic_temperature(absorbed: np.ndarray, gamma: float, config: ThermoConfig,
                          max_iter: int = 300, tolerance: float = .08) -> tuple[np.ndarray, float]:
    """Nonlinear periodic semi-infinite 1-D conduction solution in Fourier space."""
    q = np.asarray(absorbed, float)
    eps_sigma = config.emissivity*SIGMA
    # Exact instantaneous radiative equilibrium when conductivity is zero.
    if gamma == 0:
        return (np.maximum(q, 0)/eps_sigma)**.25, 0.0
    mean_t = (np.maximum(q.mean(axis=1), eps_sigma*35**4)/eps_sigma)**.25
    instant = (np.maximum(q, eps_sigma*35**4)/eps_sigma)**.25
    temperature = .65*mean_t[:, None]+.35*instant
    omega = 2*np.pi*np.fft.rfftfreq(q.shape[1], d=config.period_s/q.shape[1])
    admittance = gamma*np.sqrt(1j*omega)
    residual_max = np.inf
    for _ in range(max_iter):
        conduction = np.fft.irfft(np.fft.rfft(temperature, axis=1)*admittance[None, :],
                                  n=q.shape[1], axis=1)
        residual = q-eps_sigma*temperature**4-conduction
        residual_max = float(np.max(np.abs(residual)))
        if residual_max < tolerance:
            break
        hrad = 4*eps_sigma*np.maximum(np.mean(temperature**3, axis=1), 40**3)
        correction = np.fft.irfft(np.fft.rfft(residual, axis=1)/
                                  (hrad[:, None]+admittance[None, :]),
                                  n=q.shape[1], axis=1)
        temperature = np.clip(temperature+.42*correction, 20, 650)
    return temperature, residual_max


def _directions(config: ThermoConfig):
    phase = np.arange(config.n_phase)/config.n_phase
    angle = 2*np.pi*phase
    lat = np.deg2rad(config.subsolar_latitude_deg)
    sun = np.column_stack((np.cos(lat)*np.cos(angle), -np.cos(lat)*np.sin(angle),
                           np.full_like(angle, np.sin(lat))))
    obs_angle = angle+np.deg2rad(config.phase_angle_deg)
    obs = np.column_stack((np.cos(lat)*np.cos(obs_angle), -np.cos(lat)*np.sin(obs_angle),
                           np.full_like(angle, np.sin(lat))))
    return phase, sun, obs


def _facet_bases(normals: np.ndarray):
    ref = np.tile([0., 0., 1.], (len(normals), 1))
    ref[np.abs(normals[:, 2]) > .88] = [1., 0., 0.]
    t1 = np.cross(ref, normals)
    t1 /= np.linalg.norm(t1, axis=1)[:, None]
    t2 = np.cross(normals, t1)
    return t1, t2


def _ovirs_aperture(mesh: Mesh, observer_directions: np.ndarray,
                    observer_distances_km: np.ndarray | None,
                    boresight_directions: np.ndarray | None,
                    fov_half_angles_rad: np.ndarray | None):
    """Facet aperture weights for the public OVIRS radial-response model.

    The response is flat to 1.8 mrad and tapered to the IK radius (~2 mrad).
    The measured 4% encircled-energy wings outside the IK cylinder are included
    in the normalization but contribute space background for these observations.
    """
    if (observer_distances_km is None or boresight_directions is None or
            fov_half_angles_rad is None):
        return None
    phase_count = len(observer_directions)
    response_over_range2 = np.zeros((len(mesh.faces), phase_count))
    fill = np.zeros(phase_count)
    for phase in range(phase_count):
        spacecraft = observer_directions[phase]*observer_distances_km[phase]*1000.0
        from_spacecraft = mesh.centers-spacecraft[None, :]
        distance = np.linalg.norm(from_spacecraft, axis=1)
        look = from_spacecraft/distance[:, None]
        angle = np.arccos(np.clip(look@boresight_directions[phase], -1, 1))
        outer = float(fov_half_angles_rad[phase])
        flat = min(1.8e-3, .9*outer)
        response = np.where(angle <= flat, 1.0,
                            np.clip((outer-angle)/max(outer-flat, 1e-12), 0, 1))
        # Small-angle integral of the radial profile, enlarged for the measured
        # 4% response outside the nominal 4-mrad cylinder.
        taper_integral = ((outer**3/6 - outer*flat**2/2 + flat**3/3) /
                          max(outer-flat, 1e-12))
        omega_effective = (np.pi*flat**2+2*np.pi*taper_integral)/.96
        response_over_range2[:, phase] = response/distance**2/omega_effective
        visible_projected_solid_angle = np.sum(
            mesh.areas*np.maximum(np.einsum("ij,ij->i", mesh.normals, -look), 0)*
            response/distance**2)
        fill[phase] = visible_projected_solid_angle/omega_effective
    return response_over_range2, np.clip(fill, 0, 1)


def _crater_local_geometry(config: ThermoConfig):
    """Equal-area collocation model of a unit hemispherical crater.

    Normals point into the crater/sky and collocation positions lie on the
    spherical wall.  For a sphere, cos(theta_i)cos(theta_j)/d_ij^2 is
    1/(4R^2), giving the Rozitis view factor weights[j]/4.
    """
    cos_theta = (np.arange(config.crater_theta_bins)+.5)/config.crater_theta_bins
    theta = np.arccos(cos_theta)
    azimuth = 2*np.pi*(np.arange(config.crater_azimuth_bins)+.5)/config.crater_azimuth_bins
    local_normals = []
    for th in theta:
        for az in azimuth:
            local_normals.append([np.sin(th)*np.cos(az), np.sin(th)*np.sin(az),
                                  np.cos(th)])
    local_normals = np.asarray(local_normals, float)
    positions = -local_normals
    weights = np.full(len(local_normals), 2/len(local_normals))
    view = np.tile(weights[None, :]/4.0, (len(weights), 1))
    np.fill_diagonal(view, 0.0)
    return positions, local_normals, weights, view


def _directions_in_facet_frame(directions: np.ndarray, normals: np.ndarray,
                               t1: np.ndarray, t2: np.ndarray) -> np.ndarray:
    return np.stack((np.einsum("fc,tc->ft", t1, directions),
                     np.einsum("fc,tc->ft", t2, directions),
                     np.einsum("fc,tc->ft", normals, directions)), axis=2)


def _crater_aperture_visibility(positions: np.ndarray,
                                directions_local: np.ndarray) -> np.ndarray:
    """Ray test for whether a wall point sees a direction through the rim."""
    dz = directions_local[:, None, :, 2]
    distance = np.divide(-positions[None, :, None, 2], dz,
                         out=np.zeros((len(directions_local), len(positions),
                                       directions_local.shape[1])),
                         where=dz > 1e-12)
    exit_x = positions[None, :, None, 0] + distance*directions_local[:, None, :, 0]
    exit_y = positions[None, :, None, 1] + distance*directions_local[:, None, :, 1]
    return ((dz > 1e-12) & (distance >= 0) &
            (exit_x**2 + exit_y**2 <= 1.0 + 1e-10))


def _multiple_scattered_solar(direct_flux: np.ndarray, view: np.ndarray,
                              bond_albedo: float) -> np.ndarray:
    """Absorbed solar flux after the converged Gauss-Seidel-equivalent solve."""
    inverse = np.linalg.inv(np.eye(len(view))-bond_albedo*view)
    total_incident = np.einsum("ij,fjt->fit", inverse, direct_flux)
    return (1-bond_albedo)*total_incident


def _crater_temperatures(solar_absorbed: np.ndarray, view: np.ndarray,
                         gamma: float, config: ThermoConfig):
    shape = solar_absorbed.shape
    temperature, residual = _periodic_temperature(
        solar_absorbed.reshape(-1, shape[2]), gamma, config)
    temperature = temperature.reshape(shape)
    max_delta = np.inf
    for _ in range(max(0, config.crater_self_heating_iterations)):
        incident_thermal = (config.emissivity*SIGMA*(1-config.thermal_albedo) *
                            np.einsum("ij,fjt->fit", view, temperature**4))
        updated, residual = _periodic_temperature(
            (solar_absorbed+incident_thermal).reshape(-1, shape[2]), gamma, config)
        updated = updated.reshape(shape)
        max_delta = float(np.max(np.abs(updated-temperature)))
        temperature = 0.5*temperature+0.5*updated
        if max_delta < config.crater_self_heating_tolerance_k:
            break
    # One consistent final solve using the latest reabsorbed thermal field.
    incident_thermal = (config.emissivity*SIGMA*(1-config.thermal_albedo) *
                        np.einsum("ij,fjt->fit", view, temperature**4))
    temperature, residual = _periodic_temperature(
        (solar_absorbed+incident_thermal).reshape(-1, shape[2]), gamma, config)
    return temperature.reshape(shape), float(residual)


def simulate(mesh: Mesh, gamma: float, config: ThermoConfig = ThermoConfig(),
             sun_directions: np.ndarray | None = None,
             observer_directions: np.ndarray | None = None,
             heliocentric_distances_au: np.ndarray | None = None,
             observer_distances_km: np.ndarray | None = None,
             boresight_directions: np.ndarray | None = None,
             fov_half_angles_rad: np.ndarray | None = None) -> ModelCurves:
    """Smooth and fully cratered 4/14 um disk-radiance curves."""
    if gamma < 0:
        raise ValueError("Thermal inertia must be non-negative")
    phase, default_sun, default_obs = _directions(config)
    sun = default_sun if sun_directions is None else np.asarray(sun_directions, float)
    obs = default_obs if observer_directions is None else np.asarray(observer_directions, float)
    if sun.shape != (config.n_phase, 3) or obs.shape != (config.n_phase, 3):
        raise ValueError("SPICE direction arrays must have shape (config.n_phase, 3)")
    sun_norm = np.linalg.norm(sun, axis=1)
    obs_norm = np.linalg.norm(obs, axis=1)
    if np.any(sun_norm == 0) or np.any(obs_norm == 0):
        raise ValueError("SPICE direction vectors must be non-zero")
    sun = sun/sun_norm[:, None]
    obs = obs/obs_norm[:, None]
    if heliocentric_distances_au is None:
        solar_flux = np.full(config.n_phase,
                             SOLAR_CONSTANT/config.heliocentric_distance_au**2)
    else:
        distance = np.asarray(heliocentric_distances_au, float)
        if distance.shape != (config.n_phase,) or np.any(distance <= 0):
            raise ValueError("heliocentric_distances_au must have shape (config.n_phase,)")
        solar_flux = SOLAR_CONSTANT/distance**2
    mu_sun = np.maximum(mesh.normals@sun.T, 0)
    mu_obs = np.maximum(mesh.normals@obs.T, 0)
    if config.global_shadowing:
        from shape_radiation import directional_visibility, observer_visibility
        mu_sun *= directional_visibility(mesh, sun, mu_sun > 0)
        mu_obs *= observer_visibility(mesh, obs, observer_distances_km, mu_obs > 0)
    smooth_absorbed = (1-config.bond_albedo)*solar_flux[None, :]*mu_sun
    smooth_t, res1 = _periodic_temperature(smooth_absorbed, gamma, config)
    global_incident = np.zeros_like(smooth_absorbed)
    if config.global_self_heating:
        from shape_radiation import build_view_factors, incident_from_emitters
        edges = build_view_factors(mesh, config.global_view_factor_cache)
        for _ in range(max(1, config.global_self_heating_iterations)):
            global_incident = incident_from_emitters(
                edges, config.emissivity*SIGMA*(1-config.thermal_albedo)*smooth_t**4)
            updated, res_global = _periodic_temperature(
                smooth_absorbed+global_incident, gamma, config)
            smooth_t = .5*smooth_t+.5*updated
            res1 = max(res1, res_global)
    disk_area = np.maximum(np.sum(mesh.areas[:, None]*mu_obs, axis=0), 1e-12)

    def integrate_smooth(wavelength):
        spectral = (planck_ovirs_response(smooth_t, config.ovirs_wavelength_um,
                                          config.ovirs_fwhm_um)
                    if wavelength == 4. else planck_lambda(wavelength, smooth_t))
        return (np.sum(mesh.areas[:, None]*mu_obs*config.emissivity*
                       spectral, axis=0)/disk_area)

    s4, s14 = integrate_smooth(4.), integrate_smooth(14.)
    aperture = _ovirs_aperture(mesh, obs, observer_distances_km,
                               boresight_directions, fov_half_angles_rad)
    smooth4_detector = None
    crater4_detector_numerator = None
    modeled_fill = None
    if aperture is not None:
        aperture_weight, modeled_fill = aperture
        smooth4_detector = np.sum(
            mesh.areas[:, None]*mu_obs*config.emissivity*
            planck_ovirs_response(smooth_t, config.ovirs_wavelength_um,
                                  config.ovirs_fwhm_um)*aperture_weight, axis=0)
        crater4_detector_numerator = np.zeros(config.n_phase)
    rough4_numerator = np.zeros(config.n_phase)
    rough14_numerator = np.zeros(config.n_phase)
    res2 = 0.0
    positions, local_normals, weights, view = _crater_local_geometry(config)
    # Process the true 12k-facet mission shape without materializing all
    # facet x microfacet x phase arrays at once.
    chunk_size = max(1, int(config.facet_chunk_size))
    for first in range(0, len(mesh.faces), chunk_size):
        last = min(first + chunk_size, len(mesh.faces))
        sl = slice(first, last)
        t1, t2 = _facet_bases(mesh.normals[sl])
        sun_local = _directions_in_facet_frame(sun, mesh.normals[sl], t1, t2)
        obs_local = _directions_in_facet_frame(obs, mesh.normals[sl], t1, t2)
        sun_visible = _crater_aperture_visibility(positions, sun_local)
        raw_sun = np.maximum(np.einsum("mc,ftc->fmt", local_normals, sun_local), 0)
        raw_sun *= sun_visible
        denom_sun = np.einsum("m,fmt->ft", weights, raw_sun)
        sun_factor = np.divide(mu_sun[sl], denom_sun, out=np.zeros_like(mu_sun[sl]),
                               where=denom_sun > 0)
        direct_flux = solar_flux[None, None, :]*raw_sun*sun_factor[:, None, :]
        micro_q = _multiple_scattered_solar(direct_flux, view, config.bond_albedo)
        micro_q += global_incident[sl, None, :]
        crater_t, chunk_residual = _crater_temperatures(
            micro_q, view, gamma, config)
        res2 = max(res2, chunk_residual)
        obs_visible = _crater_aperture_visibility(positions, obs_local)
        projection = np.maximum(
            np.einsum("mc,ftc->fmt", local_normals, obs_local), 0)*obs_visible
        common = mesh.areas[sl, None, None] * weights[None, :, None] * projection * config.emissivity
        rough4_numerator += np.sum(
            common*planck_ovirs_response(crater_t, config.ovirs_wavelength_um,
                                         config.ovirs_fwhm_um), axis=(0, 1))
        if crater4_detector_numerator is not None:
            crater4_detector_numerator += np.sum(
                common*planck_ovirs_response(crater_t, config.ovirs_wavelength_um,
                                             config.ovirs_fwhm_um)*
                aperture_weight[sl, None, :], axis=(0, 1))
        rough14_numerator += np.sum(common * planck_lambda(14., crater_t), axis=(0, 1))
    r4, r14 = rough4_numerator / disk_area, rough14_numerator / disk_area
    return ModelCurves(phase, s4, r4, s14, r14, disk_area, max(res1, res2),
                       smooth4_detector, crater4_detector_numerator, modeled_fill)


def crater_fraction_to_rms_deg(fraction: float) -> float:
    return float((43/np.sqrt(.77))*np.sqrt(np.clip(fraction, 0, 1)))


def rms_deg_to_crater_fraction(rms_deg: float) -> float:
    return float(np.clip((rms_deg/(43/np.sqrt(.77)))**2, 0, 1))


def generate_grid(mesh: Mesh, gammas: Iterable[float], config: ThermoConfig):
    return {float(g): simulate(mesh, float(g), config) for g in gammas}
