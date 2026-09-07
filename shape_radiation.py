"""基于 Embree 的全局形状遮挡与宏观面元辐射交换。

本文件处理标准半球坑之外的“global shape”几何：太阳光遮挡、有限距离观测者
可见性，以及非凸宏观三角面元之间的稀疏视因子。坑内微面元过程位于
``bennu_atpm.py``，两种尺度不能混为一个视因子矩阵。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class ViewFactorEdges:
    """稀疏视因子边表。

    ``factor[k]`` 表示发射面元 ``emitter[k]`` 到接收面元
    ``receiver[k]`` 的漫射视因子；未被其它面元占据的份额记录为天空视域。
    """
    receiver: np.ndarray
    emitter: np.ndarray
    factor: np.ndarray
    sky_fraction: np.ndarray


def _intersector(mesh):
    """从项目 Mesh 构造不自动修复拓扑的 Embree 射线求交器。"""
    try:
        import trimesh
        from trimesh.ray.ray_pyembree import RayMeshIntersector
    except ImportError as exc:
        raise RuntimeError("Global ray tracing requires: pip install trimesh embreex") from exc
    tri = trimesh.Trimesh(vertices=mesh.vertices, faces=mesh.faces,
                          process=False, validate=False)
    return RayMeshIntersector(tri)


def directional_visibility(mesh, directions: np.ndarray,
                           active: np.ndarray | None = None) -> np.ndarray:
    """沿给定平行方向发射射线，返回每个面元/相位是否无遮挡。"""
    directions = np.asarray(directions, float)
    count, phases = len(mesh.faces), len(directions)
    result = np.zeros((count, phases), dtype=bool)
    intersector = _intersector(mesh)
    scale = max(float(np.max(np.linalg.norm(mesh.vertices, axis=1))), 1.0)
    # Embree 单精度遍历需要把射线起点略微推出源三角形，否则会立即击中自身。
    # 天体尺度的 1e-6 对约 500 m 的 Bennu 小于毫米，对几何结果可忽略。
    offset = scale*1e-6
    for phase, direction in enumerate(directions):
        candidates = (mesh.normals@direction > 0) if active is None else active[:, phase]
        index = np.flatnonzero(candidates)
        if not len(index):
            continue
        origins = mesh.centers[index]+mesh.normals[index]*offset
        rays = np.repeat(direction[None, :], len(index), axis=0)
        result[index, phase] = ~intersector.intersects_any(origins, rays)
    return result


def observer_visibility(mesh, observer_directions: np.ndarray,
                        observer_distances_km: np.ndarray | None,
                        active: np.ndarray | None = None) -> np.ndarray:
    """判断到有限距离航天器的可见性；距离缺失时退化为平行射线。"""
    if observer_distances_km is None:
        return directional_visibility(mesh, observer_directions, active)
    intersector = _intersector(mesh)
    count, phases = len(mesh.faces), len(observer_directions)
    result = np.zeros((count, phases), dtype=bool)
    scale = max(float(np.max(np.linalg.norm(mesh.vertices, axis=1))), 1.0)
    offset = scale*1e-6
    for phase in range(phases):
        spacecraft = observer_directions[phase]*observer_distances_km[phase]*1000.0
        delta = spacecraft[None, :]-mesh.centers
        distance = np.linalg.norm(delta, axis=1)
        rays = delta/distance[:, None]
        candidates = ((np.einsum("ij,ij->i", mesh.normals, rays) > 0) if active is None
                      else active[:, phase])
        index = np.flatnonzero(candidates)
        if not len(index):
            continue
        origins = mesh.centers[index]+mesh.normals[index]*offset
        result[index, phase] = ~intersector.intersects_any(origins, rays[index])
    return result


def build_view_factors(mesh, cache: str | Path | None = None,
                       minimum_factor: float = 1e-7,
                       chunk_size: int = 32) -> ViewFactorEdges:
    """构建含射线遮挡的宏观面元间稀疏漫射视因子。

    几何估计对应论文式(15)：``F_ij = cosθ_i cosθ_j A_j/(πd²)``。
    只有互相朝向且估计值超过阈值的面元对才做昂贵的 Embree 可见性测试。
    """
    if cache is not None and Path(cache).exists():
        data = np.load(cache)
        cache_offset = (float(data["ray_offset_fraction"])
                        if "ray_offset_fraction" in data.files else np.nan)
        if (len(data["sky_fraction"]) == len(mesh.faces) and
                np.isclose(cache_offset, 1e-6)):
            return ViewFactorEdges(data["receiver"], data["emitter"], data["factor"],
                                   data["sky_fraction"])
    intersector = _intersector(mesh)
    receivers, emitters, factors = [], [], []
    n = len(mesh.faces)
    scale = max(float(np.max(np.linalg.norm(mesh.vertices, axis=1))), 1.0)
    offset = scale*1e-6
    # 分块仅控制内存；所有满足阈值的宏观面元对仍会被检查。
    for first in range(0, n, chunk_size):
        last = min(first+chunk_size, n)
        delta = mesh.centers[None, :, :]-mesh.centers[first:last, None, :]
        distance2 = np.einsum("ijk,ijk->ij", delta, delta)
        inv_distance = np.divide(1.0, np.sqrt(distance2), out=np.zeros_like(distance2),
                                 where=distance2 > 0)
        direction = delta*inv_distance[:, :, None]
        cos_receiver = np.einsum("ik,ijk->ij", mesh.normals[first:last], direction)
        cos_emitter = -np.einsum("jk,ijk->ij", mesh.normals, direction)
        # 从接收面元 i 看发射面元 j 的视因子，面积必须取发射面元 A_j。
        estimate = (cos_receiver*cos_emitter*mesh.areas[None, :]/
                    (np.pi*np.maximum(distance2, 1e-20)))
        local_i, j = np.nonzero((cos_receiver > 0) & (cos_emitter > 0) &
                                (estimate >= minimum_factor))
        if not len(j):
            continue
        i = local_i+first
        ray_origin = mesh.centers[i]+mesh.normals[i]*offset
        ray_direction = direction[local_i, j]
        first_hit = intersector.intersects_first(ray_origin, ray_direction)
        visible = first_hit == j
        receivers.append(i[visible])
        emitters.append(j[visible])
        factors.append(estimate[local_i[visible], j[visible]])
    receiver = np.concatenate(receivers).astype(np.int32) if receivers else np.array([], np.int32)
    emitter = np.concatenate(emitters).astype(np.int32) if emitters else np.array([], np.int32)
    factor = np.concatenate(factors) if factors else np.array([], float)
    row_sum = np.bincount(receiver, weights=factor, minlength=n)
    # 面心配点近似偶尔会让某行之和略超 1。保持相对权重并为表面面元至少留下
    # 1% 天空视域，从而维持能量矩阵稳定；这是数值保护而不是新的物理参数。
    scale_rows = np.ones(n)
    excessive = row_sum > .99
    scale_rows[excessive] = .99/row_sum[excessive]
    factor *= scale_rows[receiver]
    row_sum = np.bincount(receiver, weights=factor, minlength=n)
    result = ViewFactorEdges(receiver, emitter, factor, 1-row_sum)
    if cache is not None:
        target = Path(cache)
        target.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(target, receiver=receiver, emitter=emitter, factor=factor,
                            sky_fraction=result.sky_fraction, n_faces=n,
                            ray_offset_fraction=1e-6)
    return result


def incident_from_emitters(edges: ViewFactorEdges, emitted: np.ndarray) -> np.ndarray:
    """按稀疏 ``F_ij`` 将每个发射面元/相位的出射度累加到接收面元。"""
    output = np.zeros_like(emitted)
    np.add.at(output, edges.receiver, edges.factor[:, None]*emitted[edges.emitter])
    return output
