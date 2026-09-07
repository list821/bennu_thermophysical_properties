"""Bennu 的可审计 ATPM（Advanced Thermophysical Model）正演核心。

本文件依次实现：三角形宏观面元几何、周期稳态一维导热、90° 半球坑粗糙度、
坑口遮挡、坑壁视因子自加热、太阳光多次散射、Planck 辐射及定向盘面积分。
实现依据 Rozitis & Green (2011) 公开公式独立编写，并不是作者未公开 C++
程序的逐行移植。数组下标约定通常为 ``面元 × 转相位``；粗糙坑数组则为
``宏观面元 × 坑壁微面元 × 转相位``。
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
    """不随一次正演改变的热物理、几何和数值参数。"""
    period_hours: float = 4.296061
    bond_albedo: float = 0.016
    emissivity: float = 0.9
    heliocentric_distance_au: float = 1.20
    density: float = 2500.0
    heat_capacity: float = 800.0
    phase_angle_deg: float = 5.0
    subsolar_latitude_deg: float = 0.0
    n_phase: int = 96
    # 正式计算采用 8×16=128 个等面积坑壁积分单元；早期 2×4 网格只适合调试，
    # 无法充分收敛坑口阴影边缘及定向热辐射（thermal beaming）。
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
    """宏观三角形形状模型；坐标和面积分别采用 m 与 m²。"""
    vertices: np.ndarray
    faces: np.ndarray
    centers: np.ndarray
    normals: np.ndarray
    areas: np.ndarray
    source: str


@dataclass
class ModelCurves:
    """一次热惯量正演得到的光变曲线及诊断量。"""
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
        """按论文式(24)线性混合光滑面和全坑面结果。"""
        f = float(np.clip(crater_fraction, 0.0, 1.0))
        # f 是 roughness fraction。crater_* 已用 ACF 归一到一个宏观面元的
        # 投影面积，因此这里不能再乘一次 ACF，否则会重复归一化。
        return ((1-f)*self.smooth_4um + f*self.crater_4um,
                (1-f)*self.smooth_14um + f*self.crater_14um)

    def mixed_detector(self, crater_fraction: float) -> np.ndarray:
        """式(24)在 OVIRS 孔径积分口径下的对应实现。"""
        if self.smooth_4um_detector is None or self.crater_4um_detector is None:
            return self.mixed(crater_fraction)[0]
        f = float(np.clip(crater_fraction, 0.0, 1.0))
        return (1-f)*self.smooth_4um_detector+f*self.crater_4um_detector


@dataclass(frozen=True)
class CraterGeometry:
    """单位 90° 半球坑的等面积求积几何。

    ``raw_patch_areas`` 是真实球面面积；``area_conversion_factor``（ACF）
    把坑壁总投影面积归一为一个单位宏观面元；二者的乘积
    ``area_weights`` 才是辐射积分中使用的无量纲权重。
    """

    positions: np.ndarray
    normals: np.ndarray
    raw_patch_areas: np.ndarray
    projected_area: float
    area_conversion_factor: float
    area_weights: np.ndarray
    view_factors: np.ndarray


def _mesh_from_arrays(vertices: np.ndarray, faces: np.ndarray, source: str) -> Mesh:
    """从顶点/面索引构造面心、外法向和面积，并纠正朝内的绕序。"""
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
    """读取 OBJ，多边形扇形剖分为三角形，并将 PDS 的 km 转换为 m。"""
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
    # PDS Bennu SPC OBJ 产品明确使用 km，而本地 DSK 转换文件和代理形状使用 m。
    # 这个保守尺度判断只改变长度单位，不改变法向和照明方向。
    scale_note = ""
    if np.max(np.linalg.norm(vertex_array, axis=1)) < 10.0:
        vertex_array *= 1000.0
        scale_note = " (PDS km coordinates converted to m)"
    return _mesh_from_arrays(vertex_array, np.asarray(faces, int),
                             str(obj_path.resolve()) + scale_note)


def make_bennu_proxy_mesh(n_latitude: int = 8, n_longitude: int = 16) -> Mesh:
    """生成可重复的约 246 m 半径陀螺形代理；它明确不是 SPC v13。"""
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
    """计算波长形式 Planck 辐亮度，输出单位 W m⁻² sr⁻¹ µm⁻¹。"""
    lam = np.asarray(wavelength_um, float)*1e-6
    temperature = np.maximum(np.asarray(temperature_k, float), 1.0)
    x = H*C/(lam*KB*temperature)
    return (2*H*C**2/lam**5/np.expm1(np.minimum(x, 700)))*1e-6


def planck_ovirs_response(temperature_k: np.ndarray, center_um: float,
                          fwhm_um: float, quadrature_order: int = 9) -> np.ndarray:
    """用 Gauss–Hermite 求积卷积单个 OVIRS 光谱通道的高斯线形。"""
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
    """由热惯量、体积热容和自转角频率计算昼夜热趋肤深度。"""
    return gamma/(config.density*config.heat_capacity)*np.sqrt(2/(2*np.pi/config.period_s))


def _periodic_temperature(absorbed: np.ndarray, gamma: float, config: ThermoConfig,
                          max_iter: int = 300, tolerance: float = .08) -> tuple[np.ndarray, float]:
    """在 Fourier 域求半无限介质的一维周期稳态非线性导热。

    这与论文在深度/时间网格上迭代至周期稳态所解的物理边值问题相同；差别仅是
    数值离散方式。零频项没有导热通量，非零频项的表面导热导纳为
    ``Γ sqrt(iω)``。返回表面温度和能量边界最大残差。
    """
    q = np.asarray(absorbed, float)
    eps_sigma = config.emissivity*SIGMA
    # Γ=0 时各时刻互不传热，表面温度可由瞬时辐射平衡直接求得。
    if gamma == 0:
        return (np.maximum(q, 0)/eps_sigma)**.25, 0.0
    mean_t = (np.maximum(q.mean(axis=1), eps_sigma*35**4)/eps_sigma)**.25
    instant = (np.maximum(q, eps_sigma*35**4)/eps_sigma)**.25
    temperature = .65*mean_t[:, None]+.35*instant
    omega = 2*np.pi*np.fft.rfftfreq(q.shape[1], d=config.period_s/q.shape[1])
    admittance = gamma*np.sqrt(1j*omega)
    residual_max = np.inf
    for _ in range(max_iter):
        # 向地下为正的传导通量；能量边界为 q_abs = εσT⁴ + q_cond。
        conduction = np.fft.irfft(np.fft.rfft(temperature, axis=1)*admittance[None, :],
                                  n=q.shape[1], axis=1)
        residual = q-eps_sigma*temperature**4-conduction
        residual_max = float(np.max(np.abs(residual)))
        if residual_max < tolerance:
            break
        # 用辐射项的一阶导数作预条件，再以 0.42 欠松弛避免非线性振荡。
        hrad = 4*eps_sigma*np.maximum(np.mean(temperature**3, axis=1), 40**3)
        correction = np.fft.irfft(np.fft.rfft(residual, axis=1)/
                                  (hrad[:, None]+admittance[None, :]),
                                  n=q.shape[1], axis=1)
        temperature = np.clip(temperature+.42*correction, 20, 650)
    return temperature, residual_max


def _directions(config: ThermoConfig):
    """没有逐帧 SPICE 输入时，生成等间隔的太阳与观测者单位方向。"""
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
    """为每个宏观面元建立右手正交局部基 ``(t1,t2,n)``。"""
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
    """依据公开 OVIRS 径向响应计算各宏观面元的孔径权重。

    响应在 1.8 mrad 内取常数，随后衰减至 IK 给出的约 2 mrad 半径。IK 圆柱
    外约 4% 的包围能量计入归一化；对本次观测而言该部分落在太空背景上。
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
        # 小角近似下积分径向响应，再为名义 4 mrad 圆柱外的 4% 能量修正。
        taper_integral = ((outer**3/6 - outer*flat**2/2 + flat**3/3) /
                          max(outer-flat, 1e-12))
        omega_effective = (np.pi*flat**2+2*np.pi*taper_integral)/.96
        response_over_range2[:, phase] = response/distance**2/omega_effective
        visible_projected_solid_angle = np.sum(
            mesh.areas*np.maximum(np.einsum("ij,ij->i", mesh.normals, -look), 0)*
            response/distance**2)
        fill[phase] = visible_projected_solid_angle/omega_effective
    return response_over_range2, np.clip(fill, 0, 1)


def _hemispherical_crater_geometry(config: ThermoConfig) -> CraterGeometry:
    """生成论文所用 90° spherical-section crater 的等面积求积网格。

    半球半径取 1、坑口位于 z=0。法向指向坑腔/天空，配点位于球壁。这里没有
    改成显式三角片，是因为等面积中点求积已经覆盖同一连续半球；128 个积分
    单元仅在阴影边界处存在有限分辨率误差。

    ACF 按论文附录定义为：单位宏观面元面积除以全部坑壁单元在局部水平面上的
    总投影面积。对单位半球，该投影面积应为 π，故 ACF=1/π，最终每个等面积
    单元的式(24)权重为 ``(1/π)(2π/N)=2/N``。把这一关系显式保存可防止后续
    再乘 ACF 而造成重复归一化。
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
    count = len(local_normals)
    raw_patch_areas = np.full(count, 2*np.pi/count)
    projected_area = float(np.sum(raw_patch_areas*local_normals[:, 2]))
    area_conversion_factor = 1.0/projected_area
    weights = area_conversion_factor*raw_patch_areas

    # 球面任意两个不同微面元间，cosθ_i cosθ_j/d² 恒为 1/(4R²)。再计入
    # 式(15)的 1/π 后，F_ij=A_j/(4πR²)。单位半球的 ACF=1/π，故恰好写成
    # ACF×A_j/4 = weights[j]/4。
    view = np.tile(weights[None, :]/4.0, (count, 1))
    np.fill_diagonal(view, 0.0)
    return CraterGeometry(positions, local_normals, raw_patch_areas,
                          projected_area, area_conversion_factor, weights, view)


def _crater_local_geometry(config: ThermoConfig):
    """兼容旧调用的四元组接口；数值与显式 ACF 重构前完全一致。"""
    crater = _hemispherical_crater_geometry(config)
    return (crater.positions, crater.normals, crater.area_weights,
            crater.view_factors)


def _directions_in_facet_frame(directions: np.ndarray, normals: np.ndarray,
                               t1: np.ndarray, t2: np.ndarray) -> np.ndarray:
    """把惯性/天体固定系方向投影到每个宏观面元的局部坐标。"""
    return np.stack((np.einsum("fc,tc->ft", t1, directions),
                     np.einsum("fc,tc->ft", t2, directions),
                     np.einsum("fc,tc->ft", normals, directions)), axis=2)


def _crater_aperture_visibility(positions: np.ndarray,
                                directions_local: np.ndarray) -> np.ndarray:
    """由射线与 z=0 坑口平面的交点判断坑壁微面元是否看见目标方向。"""
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
    """求解 ``(I-AF)⁻¹``，得到坑内多次散射后的吸收太阳通量。"""
    inverse = np.linalg.inv(np.eye(len(view))-bond_albedo*view)
    total_incident = np.einsum("ij,fjt->fit", inverse, direct_flux)
    return (1-bond_albedo)*total_incident


def _crater_temperatures(solar_absorbed: np.ndarray, view: np.ndarray,
                         gamma: float, config: ThermoConfig):
    """迭代坑壁热辐射自加热与一维导热，返回各微面元周期温度。"""
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
    # 用最后一次再吸收热辐射场作一致的终解，避免返回混合迭代的中间温度。
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
    """计算同一形状在全光滑与全半球坑两种端元下的 4/14 µm 光变。

    外部传入的形状、SPICE 几何、距离和 OVIRS 孔径输入原样使用。本函数只做
    正演，不估计粗糙度比例；反演阶段再通过 :meth:`ModelCurves.mixed` 按
    论文式(24)组合两个端元。
    """
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
    # 逐帧日心距优先；缺省时才使用配置中的单一平均值。
    if heliocentric_distances_au is None:
        solar_flux = np.full(config.n_phase,
                             SOLAR_CONSTANT/config.heliocentric_distance_au**2)
    else:
        distance = np.asarray(heliocentric_distances_au, float)
        if distance.shape != (config.n_phase,) or np.any(distance <= 0):
            raise ValueError("heliocentric_distances_au must have shape (config.n_phase,)")
        solar_flux = SOLAR_CONSTANT/distance**2
    # 宏观面元的太阳/观测余弦；背光或背向观测者时截为 0。
    mu_sun = np.maximum(mesh.normals@sun.T, 0)
    mu_obs = np.maximum(mesh.normals@obs.T, 0)
    if config.global_shadowing:
        # 对非凸全局形状追加 Embree 射线遮挡；不改变入射/观测几何本身。
        from shape_radiation import directional_visibility, observer_visibility
        mu_sun *= directional_visibility(mesh, sun, mu_sun > 0)
        mu_obs *= observer_visibility(mesh, obs, observer_distances_km, mu_obs > 0)
    # 光滑端元：吸收太阳通量 → 周期稳态表面温度。
    smooth_absorbed = (1-config.bond_albedo)*solar_flux[None, :]*mu_sun
    smooth_t, res1 = _periodic_temperature(smooth_absorbed, gamma, config)
    global_incident = np.zeros_like(smooth_absorbed)
    if config.global_self_heating:
        # 非凸宏观面元间的长波辐射交换。该迭代只作用于全局形状端元。
        from shape_radiation import build_view_factors, incident_from_emitters
        edges = build_view_factors(mesh, config.global_view_factor_cache)
        for _ in range(max(1, config.global_self_heating_iterations)):
            global_incident = incident_from_emitters(
                edges, config.emissivity*SIGMA*(1-config.thermal_albedo)*smooth_t**4)
            updated, res_global = _periodic_temperature(
                smooth_absorbed+global_incident, gamma, config)
            smooth_t = .5*smooth_t+.5*updated
            res1 = max(res1, res_global)
    # 可见宏观面元的投影面积，用于把盘积分功率换成盘平均辐亮度。
    disk_area = np.maximum(np.sum(mesh.areas[:, None]*mu_obs, axis=0), 1e-12)

    def integrate_smooth(wavelength):
        """对全光滑端元作可见投影面积加权的方向积分。"""
        spectral = (planck_ovirs_response(smooth_t, config.ovirs_wavelength_um,
                                          config.ovirs_fwhm_um)
                    if wavelength == 4. else planck_lambda(wavelength, smooth_t))
        return (np.sum(mesh.areas[:, None]*mu_obs*config.emissivity*
                       spectral, axis=0)/disk_area)

    s4, s14 = integrate_smooth(4.), integrate_smooth(14.)
    # OVIRS 检测器口径积分与上面的盘平均辐亮度是两种不同口径；二者都保留，
    # 防止把已是检测器辐亮度的数据再除一次 FOV 填充率。
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
    crater = _hemispherical_crater_geometry(config)
    positions = crater.positions
    local_normals = crater.normals
    weights = crater.area_weights
    view = crater.view_factors
    # 对 12k 面元任务形状分块，避免一次构造完整的“宏观面元×坑壁×相位”数组。
    chunk_size = max(1, int(config.facet_chunk_size))
    for first in range(0, len(mesh.faces), chunk_size):
        last = min(first + chunk_size, len(mesh.faces))
        sl = slice(first, last)
        # 将同一个标准半球坑旋转到每个宏观面元的局部切平面。
        t1, t2 = _facet_bases(mesh.normals[sl])
        sun_local = _directions_in_facet_frame(sun, mesh.normals[sl], t1, t2)
        obs_local = _directions_in_facet_frame(obs, mesh.normals[sl], t1, t2)
        # 坑壁太阳直射：入射余弦与坑口射线可见性共同决定直接照明。
        sun_visible = _crater_aperture_visibility(positions, sun_local)
        raw_sun = np.maximum(np.einsum("mc,ftc->fmt", local_normals, sun_local), 0)
        raw_sun *= sun_visible
        # 有限坑壁求积会使投影和宏观面元余弦略有偏差；sun_factor 强制能量守恒，
        # 使全坑端元通过坑口接收的直射总功率等于对应宏观面元。
        denom_sun = np.einsum("m,fmt->ft", weights, raw_sun)
        sun_factor = np.divide(mu_sun[sl], denom_sun, out=np.zeros_like(mu_sun[sl]),
                               where=denom_sun > 0)
        direct_flux = solar_flux[None, None, :]*raw_sun*sun_factor[:, None, :]
        # 太阳光多次散射、可能的全局自加热以及坑壁之间的长波自加热依次进入边界。
        micro_q = _multiple_scattered_solar(direct_flux, view, config.bond_albedo)
        micro_q += global_incident[sl, None, :]
        crater_t, chunk_residual = _crater_temperatures(
            micro_q, view, gamma, config)
        res2 = max(res2, chunk_residual)
        # thermal beaming：只有能从坑口看见且朝向观测者的热微面元才贡献辐射。
        obs_visible = _crater_aperture_visibility(positions, obs_local)
        projection = np.maximum(
            np.einsum("mc,ftc->fmt", local_normals, obs_local), 0)*obs_visible
        # weights = ACF × 坑壁真实面积。这里就是论文式(24)粗糙端元内层求和；
        # 因此外层 mixed() 只需乘 roughness fraction，无需再乘 ACF。
        common = (mesh.areas[sl, None, None] * weights[None, :, None] *
                  projection * config.emissivity)
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
    """按本项目采用的论文标定把半球坑覆盖率换为 RMS 坡度。"""
    return float((43/np.sqrt(.77))*np.sqrt(np.clip(fraction, 0, 1)))


def rms_deg_to_crater_fraction(rms_deg: float) -> float:
    """把 RMS 坡度换为半球坑覆盖率，并限制在物理区间 [0,1]。"""
    return float(np.clip((rms_deg/(43/np.sqrt(.77)))**2, 0, 1))


def generate_grid(mesh: Mesh, gammas: Iterable[float], config: ThermoConfig):
    """对一组热惯量逐一执行正演，供后续 χ² 网格搜索调用。"""
    return {float(g): simulate(mesh, float(g), config) for g in gammas}
