# Bennu 逐帧 SPICE 几何说明

## SPICE 是什么

SPICE 是 NASA/JPL NAIF 维护的任务辅助几何系统，名称来自 **S**pacecraft、
**P**lanet、**I**nstrument、**C**-matrix、**E**vents。它把不同核文件中的时间、
轨道、天体自转、航天器姿态、仪器视场和事件定义组合起来，在指定 UTC 时刻
计算一致的空间几何。

SPICE 不计算表面温度，也不提供 OVIRS 辐亮度。它在本项目中回答的是：
“这一帧曝光时，太阳和航天器分别位于 Bennu 哪个方向、距离多远？”这些量再
送入 ATPM 的日照、可见性、导热和盘积分模块。

## 当前使用的核

|核类型|文件|作用|
|---|---|---|
|LSK|`naif0012.tls`|UTC、闰秒和 ET/TDB 转换|
|PCK|`pck00010.tpc`, `bennu_v14.tpc`|Bennu 自转轴、本初子午线、自转率；选 v14 以配准 SPCv14|
|FK|`orx_v14.tf`, `orx_shape_v03.tf`|OSIRIS-REx 坐标系和形状表面编号|
|SPK|`de424.bsp`|太阳系天体轨道|
|SPK|`bennu_refdrmc_v1.bsp`|Bennu 轨道|
|SPK|`orx_180801_190302_181218_od077_v1.bsp`|覆盖 2018-11-02/03 的重建航天器轨道|
|SPK|`orx_struct_v04.bsp`|OSIRIS-REx 结构位置定义|
|SCLK|`orx_sclkscet_00093.tsc`|航天器时钟与 ET 的转换|
|CK|`orx_sc_rel_181029_181104_v02.bc`|覆盖两观测日的重建航天器姿态|
|IK|`orx_ovirs_v00.ti`|OVIRS 科学视轴、圆形视场和边界|

当前已加载 CK、SCLK 和 IK。IK 的科学通道 NAIF ID 为 `-64321`、仪器系为
`ORX_OVIRS_SCI`、视轴为 `+Z`，圆形半角为约 0.115°（约 2.007 mrad）。

## 每帧计算

对每个 FITS 的 `MID_OBS` UTC：

1. LSK 把 UTC 转成 SPICE ET。
2. `spkpos` 在 `IAU_BENNU` 固连系中求 Bennu→太阳和 Bennu→OSIRIS-REx 向量，
   使用 `LT+S` 光行时与恒星像差修正。
3. `pxform('ORX_OVIRS_SCI', 'IAU_BENNU', et)` 把 CK 给出的 OVIRS `+Z` 视轴
   旋转到 Bennu 固连系；`getfov(-64321)` 读取 IK 的圆形视场边缘。
4. 两位置向量单位化为 `sun_[xyz]` 和 `observer_[xyz]`，同时保存
   `boresight_[xyz]` 与有限距离的航天器位置。
5. 计算日心距离、航天器距离、相位角
   `acos(sun_direction · observer_direction)` 和太阳下点纬度。
6. ATPM 对形状面元法向 `n` 使用
   `max(0, n · sun_direction)` 求入射余弦、
   `max(0, n · observer_direction)` 求可见投影，太阳通量按 `1361/r_au²` 缩放。

全部 34,254 帧先独立求几何，再按论文要求和辐亮度一起分到 1° 自转相位箱；
热传导求解时将这些方向周期插值到统一相位网格。既然 SPICE 与 SPC 已给出绝对
经度，本流程不再额外拟合任意自转相位平移。

## 数值核验

|UTC|日心距离 (au)|航天器距离 (km)|相位角|
|---|---:|---:|---:|
|2018-11-02 04:12|1.040938|197.749|5.1806°|
|2018-11-03 04:12|1.037625|191.195|4.5768°|

这与论文的约 1.041/1.037 au、197/190 km、5.1°/4.5°相符。

首帧 CK/IK 视轴为 RA=345.3289269°、Dec=-7.4104782°；FITS 头为
RA=345.3289294°、Dec=-7.4104751°，相差约 4×10⁻⁶°，验证了 CK、SCLK 与帧时间
的配准。

## 运行

```powershell
cd D:\ATPM
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' .\download_spice_kernels.py
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -u .\reproduce_supplementary_fig2b.py
```

`D:\ATPM` 是指向实际项目的 ASCII 目录联接。Windows 版 CSPICE N0067 不能可靠
加载含中文字符的核路径，因此运行 SPICE 流程时应从该入口进入。

逐帧结果位于 `output/supplementary_fig2b/processed_ovirs_frames.csv`，1° 分箱后的
辐亮度与方向位于 `output/supplementary_fig2b/ovirs_thermal_1deg.csv`。
