"""Auditable, dependency-light thermophysical model for asteroid Bennu.

This implements triangular-facet illumination, periodic 1-D conduction,
fractional hemispherical-crater roughness, Planck radiation and disk integration.
It is not the unpublished 2019 author code: the built-in mesh is a Bennu-like
proxy and crater self-heating/ray tracing are documented approximations.
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
    crater_theta_bins: int = 3
    crater_azimuth_bins: int = 6

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
    max_boundary_residual_w_m2: float

    def mixed(self, crater_fraction: float) -> tuple[np.ndarray, np.ndarray]:
        f = float(np.clip(crater_fraction, 0.0, 1.0))
        return ((1-f)*self.smooth_4um + f*self.crater_4um,
                (1-f)*self.smooth_14um + f*self.crater_14um)


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
    return _mesh_from_arrays(np.asarray(vertices, float), np.asarray(faces, int),
                             str(obj_path.resolve()))


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


def planck_lambda(wavelength_um: float, temperature_k: np.ndarray) -> np.ndarray:
    lam = float(wavelength_um)*1e-6
    temperature = np.maximum(np.asarray(temperature_k, float), 1.0)
    x = H*C/(lam*KB*temperature)
    return (2*H*C**2/lam**5/np.expm1(np.minimum(x, 700)))*1e-6


def skin_depth_m(gamma: float, config: ThermoConfig) -> float:
    return gamma/(config.density*config.heat_capacity)*np.sqrt(2/(2*np.pi/config.period_s))


def _periodic_temperature(absorbed: np.ndarray, gamma: float, config: ThermoConfig,
                          max_iter: int = 100, tolerance: float = .08) -> tuple[np.ndarray, float]:
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


def _crater_micro_normals(normals: np.ndarray, config: ThermoConfig):
    ref = np.tile([0., 0., 1.], (len(normals), 1))
    ref[np.abs(normals[:, 2]) > .88] = [1., 0., 0.]
    t1 = np.cross(ref, normals)
    t1 /= np.linalg.norm(t1, axis=1)[:, None]
    t2 = np.cross(normals, t1)
    cos_theta = (np.arange(config.crater_theta_bins)+.5)/config.crater_theta_bins
    theta = np.arccos(cos_theta)
    azimuth = 2*np.pi*(np.arange(config.crater_azimuth_bins)+.5)/config.crater_azimuth_bins
    micro = []
    for th in theta:
        for az in azimuth:
            micro.append(np.cos(th)*normals+np.sin(th)*(np.cos(az)*t1+np.sin(az)*t2))
    result = np.stack(micro, axis=1)
    # Actual hemisphere area / aperture area = 2.
    return result, np.full(result.shape[1], 2/result.shape[1])


def simulate(mesh: Mesh, gamma: float, config: ThermoConfig = ThermoConfig()) -> ModelCurves:
    """Smooth and fully cratered 4/14 um disk-radiance curves."""
    if gamma < 0:
        raise ValueError("Thermal inertia must be non-negative")
    phase, sun, obs = _directions(config)
    scale = (1-config.bond_albedo)*SOLAR_CONSTANT/config.heliocentric_distance_au**2
    mu_sun = np.maximum(mesh.normals@sun.T, 0)
    mu_obs = np.maximum(mesh.normals@obs.T, 0)
    smooth_t, res1 = _periodic_temperature(scale*mu_sun, gamma, config)
    micro, weights = _crater_micro_normals(mesh.normals, config)
    raw_sun = np.maximum(np.einsum("fmc,tc->fmt", micro, sun), 0)
    denom_sun = np.einsum("m,fmt->ft", weights, raw_sun)
    sun_factor = np.divide(mu_sun, denom_sun, out=np.zeros_like(mu_sun), where=denom_sun>0)
    micro_q = scale*raw_sun*sun_factor[:, None, :]
    crater_flat, res2 = _periodic_temperature(micro_q.reshape(-1, config.n_phase), gamma, config)
    crater_t = crater_flat.reshape(len(mesh.faces), len(weights), config.n_phase)
    raw_obs = np.maximum(np.einsum("fmc,tc->fmt", micro, obs), 0)
    denom_obs = np.einsum("m,fmt->ft", weights, raw_obs)
    obs_factor = np.divide(mu_obs, denom_obs, out=np.zeros_like(mu_obs), where=denom_obs>0)
    projection = raw_obs*obs_factor[:, None, :]
    disk_area = np.maximum(np.sum(mesh.areas[:, None]*mu_obs, axis=0), 1e-12)

    def integrate(wavelength):
        smooth = np.sum(mesh.areas[:, None]*mu_obs*config.emissivity*
                        planck_lambda(wavelength, smooth_t), axis=0)/disk_area
        rough = np.sum(mesh.areas[:, None, None]*weights[None, :, None]*projection*
                       config.emissivity*planck_lambda(wavelength, crater_t), axis=(0, 1))/disk_area
        return smooth, rough
    s4, r4 = integrate(4.)
    s14, r14 = integrate(14.)
    return ModelCurves(phase, s4, r4, s14, r14, max(res1, res2))


def crater_fraction_to_rms_deg(fraction: float) -> float:
    return float((43/np.sqrt(.77))*np.sqrt(np.clip(fraction, 0, 1)))


def rms_deg_to_crater_fraction(rms_deg: float) -> float:
    return float(np.clip((rms_deg/(43/np.sqrt(.77)))**2, 0, 1))


def generate_grid(mesh: Mesh, gammas: Iterable[float], config: ThermoConfig):
    return {float(g): simulate(mesh, float(g), config) for g in gammas}
