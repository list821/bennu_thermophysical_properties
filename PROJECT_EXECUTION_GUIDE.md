# Bennu ATPM Python 项目：运行、输入、代码与结果说明

## 1. 项目目标与结论边界

本项目用公开 OSIRIS-REx 数据重建 Bennu 的 4 μm OVIRS 自转光变，并用 ATPM 风格的热物理模型扫描热惯量 Γ。核心链路是：

```text
OVIRS L2 FITS ──太阳反射扣除──> 每帧 4 μm 热辐亮度
       │                              │
       └──时间──> SPICE SPK/CK/IK ───┤
                                      v
SPCv14 OBJ ──面元/光线/视因子──> 逐帧几何与 1° 相位分箱
                                      │
                                      v
              一维导热 + 坑粗糙度 + 自加热 + FOV 响应
                                      │
                                      v
                         Γ 网格模型曲线与 χ²
```

当前公开输入下的审计扫描最优值为 Γ=248 J m⁻² K⁻¹ s⁻½，但约化 χ²=154.24，统计上不可接受。因此它是“当前公开数据重建的数值最小值”，不是对论文 Γ=350±20 的独立替代测量。

## 2. 环境与目录

项目实际目录为 `D:\ChatGPT\热物理反演`。Windows CSPICE N0067 对含中文的核路径兼容性不好，因此 SPICE 命令从 ASCII 联接 `D:\ATPM` 运行。

Python 依赖在 `requirements_atpm.txt`：

- `numpy`：数组、FFT、线性代数和数值积分；
- `spiceypy`：NASA SPICE 的 Python 接口；
- `trimesh`：三角网格对象；
- `embreex`：Embree 加速的射线求交。

安装依赖：

```powershell
cd D:\ATPM
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -m pip install -r .\requirements_atpm.txt
```

## 3. 输入数据逐项解释

### 3.1 OVIRS L2 v2 FITS

下载目录：`data/ovirs/full/`。下载器只选择 2018-11-02 和 2018-11-03 的 `*_ovr_scil2_calv2.fits`，共 34,254 帧，每帧 558,720 字节。

每个 FITS 的固定布局为：

|内容|数组形状|类型/字节序|文件偏移|含义|
|---|---:|---|---:|---|
|Primary header|5760 bytes|80 字节 FITS card|0|时间、填充率、距离、指向等标量元数据|
|radiance|20×512|`>f8`|5760|20 个 superpixel、每个 512 个光谱样本的标定辐亮度|
|quality|20×512|`>i4`|92160|逐样本质量位；低 4 位是有效探测器像元数，bit 4 为空，bit 6 为 outlier|
|wavelength|20×512|`>f8`|138240|每个辐亮度样本对应的波长，单位 μm|
|uncertainty|20×512|`>f8`|475200|逐样本的 1σ 辐亮度不确定度|

代码优先读取同名 PDS4 XML；若 XML 不在本地，则从 FITS 头读取：

- `MIDOBS`：曝光中点 UTC，是 SPICE 几何的时间输入；
- `BS_FLAG`：视轴/目标状态标志；
- `FILL_FAC`：OVIRS 视场与 Bennu 投影相交的几何比例；
- `PHASEANG`：太阳—Bennu—航天器相角；
- `SUN_RNG`：日心距离；
- `TARGRNG`：航天器到目标距离。

注意：`FILL_FAC` 不是温度，也不是辐亮度校准系数。它描述目标填充视场的比例。本项目同时保留：

- `thermal_radiance_detector`：直接从 L2 扣除反射光后的探测器视场平均量；
- `thermal_radiance = detector/FILL_FAC`：用于和盘面平均模型比较的盘面等效量。

### 3.2 OSIRIS-REx 官方太阳光谱

文件：`data/ovirs/orexsolarflux.csv`。列为波长 μm、太阳光谱能量密度、误差和光谱平滑分辨率。它用于给反射太阳光提供光谱形状：

1. 在 2.05–2.15 μm 选择有效样本；
2. 用不确定度倒数平方作为权重拟合尺度 `a`；
3. 在 3.98–4.02 μm 计算 `L_thermal = L_measured - a L_solar`；
4. 对剩余样本求均值和误差。

2.1 μm 是反射光锚点；4 μm 附近热辐射占主导。不能把整个 3.5–4.0 μm 的平均辐亮度直接和单色 4 μm 模型比较，因为该区间位于热谱 Wien 侧，斜率很陡。

### 3.3 SPCv14 三角形状

文件：`data/shape/g_12560mm_spc_obj_0000n00000_v014.obj`，包含 6,534 个顶点和 12,288 个三角面元。

OBJ 文本语法：

```text
v x y z
f i j k
```

`v` 是顶点坐标，PDS 产品单位为 km；`f` 是从 1 开始的顶点编号。`load_obj()` 自动把半径小于 10 的 PDS 坐标乘 1000 转成 m，并计算：

- 三角中心 `centers`；
- 外法向 `normals`；
- 面积 `areas`；
- 方向错误的三角形会交换两个顶点以翻转法向。

论文使用 SPCv13；公开可直接获得的最近替代是 SPCv14。这会产生面元地形、法向和自转经度上的小差别。

### 3.4 SPICE 核

`download_spice_kernels.py` 下载：

|核|作用|
|---|---|
|LSK `naif0012.tls`|UTC 与 ET/TDB、闰秒|
|PCK `pck00010.tpc`, `bennu_v14.tpc`|天体常数、Bennu 自转轴/本初子午线/自转率|
|FK `orx_v14.tf`, `orx_shape_v03.tf`|任务坐标系定义|
|SPK `de424.bsp`, Bennu/ORX SPK|太阳、Bennu 和航天器的位置|
|SCLK `orx_sclkscet_00093.tsc`|航天器时钟到 ET 的映射|
|CK `orx_sc_rel_181029_181104_v02.bc`|两次观测期间的重建姿态|
|IK `orx_ovirs_v00.ti`|OVIRS `+Z` 视轴、圆形视场及边界|

每帧输出：太阳方向、航天器方向和位置、OVIRS 视轴、FOV 半角、日心/航天器距离、相角和太阳下点纬度，全部表达在 `IAU_BENNU` 固连坐标系。

### 3.5 模型参数

`ThermoConfig` 的主要参数：

|字段|默认值|物理意义|
|---|---:|---|
|`period_hours`|4.296061 h|Bennu 自转周期|
|`bond_albedo`|0.016|Bond 反照率|
|`emissivity`|0.9|热红外发射率|
|`density`|2500 kg m⁻³|材料密度，仅在趋肤深度换算使用|
|`heat_capacity`|800 J kg⁻¹ K⁻¹|比热|
|`n_phase`|48–96|一个自转周期的数值相位点数|
|`crater_theta_bins × crater_azimuth_bins`|8×16|每个半球坑的 128 个微面元|
|`crater_self_heating_iterations`|4|坑内热辐射迭代次数|
|`ovirs_wavelength_um`|4.00038 μm|模型 OVIRS 通道中心|
|`ovirs_fwhm_um`|0.00641 μm|公开采样推导的高斯响应 FWHM|
|`global_shadowing`|运行时开启|全局形状太阳/观测遮挡|
|`global_self_heating`|运行时开启|全局面元相互热辐射|

## 4. 推荐的完整运行顺序

### 第一步：下载公开输入

```powershell
cd D:\ATPM
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' .\download_spice_kernels.py
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -u .\download_ovirs_20181102_03.py --workers 24 --output .\data\ovirs\full
```

下载器使用 `.part` 临时文件，完成后原子改名；已有且大小正确的 FITS 会跳过，因此可以断点续传。

### 第二步：测试

```powershell
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -m unittest -v .\test_bennu_atpm.py .\test_pds_observations.py
```

目前 7 项测试通过，覆盖普朗克函数、响应卷积、粗糙度换算、热边界残差、坑视因子/遮挡、OTES 二进制和 OVIRS FITS 读取。

### 第三步：真实 OVIRS 反演

```powershell
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  .\reproduce_supplementary_fig2b.py `
  --gamma-min 240 --gamma-max 256 --gamma-step 2 `
  --model-phases 48 `
  --fit-facet-stride 128 `
  --facet-chunk-size 16 `
  --crater-theta-bins 8 --crater-azimuth-bins 16 `
  --skip-full-shape-figure
```

`--fit-facet-stride 128` 表示拟合阶段每 128 个 SPCv14 面元保留一个，共 96 个；保留面元面积乘 128 以维持投影面积。它是速度—精度折中，不等于完整 12,288 面元反演。去掉 `--skip-full-shape-figure` 可用全部面元绘图，但当前机器六条参考曲线约需 40 分钟。

### 第四步：阅读结果

- `summary.json`：模型配置、最优 Γ、χ²、论文参考与限制；
- `ovirs_gamma_phase_chi2.csv`：每个 Γ 的 χ²；
- `ovirs_thermal_1deg.csv`：两天各 360 个相位箱的探测器/盘面辐亮度、误差、填充率和几何；
- `supplementary_figure_2b_reproduction.svg`：真实数据与 Γ=330/350/370 模型曲线。

## 5. 各 Python 文件和关键语法

### `bennu_atpm.py`：热物理正演核心

- `@dataclass(frozen=True)`：自动生成初始化器；`frozen=True` 防止配置在计算中被意外修改。
- `float | np.ndarray`：Python 类型联合，表示标量或 NumPy 数组均可。
- `np.einsum()`：爱因斯坦求和，例如 `ij,ij->i` 是逐行点积；避免显式 Python 循环。
- `[:, None]`：增加长度为 1 的轴，利用广播把逐相位量乘到全部面元。
- `np.fft.rfft/irfft`：在频域求周期稳态地下导热。

核心方程：

```text
(1-A) S/r² μ + Q_self = εσT⁴ + Q_cond
Q_cond,n = Γ sqrt(i n ω) T_n
```

`_periodic_temperature()` 对非线性辐射边界迭代；`planck_lambda()` 计算光谱辐亮度；`planck_ovirs_response()` 用 9 点 Gauss–Hermite 积分卷积高斯响应；`simulate()` 把光滑面和坑面按覆盖率混合并作盘积分/FOV 积分。

### `shape_radiation.py`：全局遮挡与自加热

- `try/except ImportError`：依赖缺失时给出明确安装错误；
- `intersects_any()`：判断太阳/观测射线是否被任何面元挡住；
- `intersects_first()`：判断接收面元发出的射线首先击中哪个发射面元；
- `np.add.at()`：根据稀疏接收面元编号累加辐射贡献。

视因子近似为：

```text
F_ij = cosθ_i cosθ_j A_j / (π r_ij²)
```

再通过 Embree 光线追踪删除不可互见的面元对。完整 SPCv14 缓存有 220,623 条可见交换边。

### `pds_observations.py`：PDS 数据解析

- `Path.open('rb')`：二进制读取；
- `'>f8'` / `'>i4'`：大端 64 位浮点/32 位整数；
- `quality & 0x40`：位掩码语法，用于检查质量标志；
- `with ... as stream`：离开代码块时自动关闭文件；
- 加权尺度使用 `w=1/σ²`，反射扣除后误差取传播误差和样本散布标准误差的较大者。

### `spice_geometry.py`：逐帧任务几何

- `try/finally` 保证即使计算失败也执行 `spice.kclear()`；
- `spkpos()` 求相对位置，`LT+S` 表示光行时和恒星像差修正；
- `getfov(-64321, 4)` 读取 OVIRS IK；
- `pxform()` 把 `ORX_OVIRS_SCI` 视轴旋转到 `IAU_BENNU`；
- 字典 `return {...}` 使每一种几何数组有明确字段名。

### `reproduce_supplementary_fig2b.py`：当前主入口

执行顺序是 `process_frames → attach_spice_geometry → bin_one_degree → sampled_mesh → fit_gamma → write_figure`。

- `argparse` 把命令行参数转成 `args.gamma_min` 等字段；
- `np.mod((t-t0)/P,1)` 把 UTC 映射为 0–1 自转相位；
- `floor(phase*360)` 得到 1° 箱编号；
- 每箱用 5×MAD 去除离群点；
- χ² 为 `Σ[(model-data)/σ]²`；SPICE 已给出绝对经度，因此不再自由平移相位。

### 下载和辅助脚本

- `download_ovirs_20181102_03.py`：并行下载并校验两天 OVIRS FITS；
- `download_spice_kernels.py`：下载最小 SPICE 核集合；
- `download_paper_inputs.py`：较早的综合输入下载器，包含 OTES 和替代形状；
- `convert_bennu_dsk_to_obj.py`：通过 SPICE 把 DSK 三角面转换为 OBJ；
- `crater_resolution_convergence.py`：比较不同坑网格分辨率；
- `reproduce_bennu.py`：用论文参数生成代理观测的自洽测试，不是真实数据反演；
- `invert_real_bennu.py`：早期 OVIRS/OTES 联合实验入口，保留用于比较；当前补充图 2b 应使用 `reproduce_supplementary_fig2b.py`；
- `test_bennu_atpm.py`, `test_pds_observations.py`：单元测试。

## 6. 当前结果和论文差异

精细网格：

|Γ|χ²|
|---:|---:|
|240|125125.72|
|244|113474.65|
|246|111034.81|
|248|**110746.50**|
|250|112534.17|
|252|116324.92|
|256|129634.62|

最优 Γ=248，但与论文 350±20 有显著差异，原因按重要性分为：

1. **作者 L3a 未公开**：公开档案是 L2 v2。校准文档明确提示目标欠填充 FOV 会产生光谱伪影，但作者使用的逐帧/逐波长经验校正数组没有公开。`FILL_FAC` 只能修正几何稀释，不能恢复该经验校正。
2. **观测/模型归一化曾混用**：探测器视场平均约 1.12×10⁻⁴，而除以填充率的盘面等效约 3.07×10⁻⁴。原先所谓“模型高两倍”主要是将前者和盘面模型比较。统一定义后，Γ=350 模型反而低约 29%。
3. **SPCv14 代替 SPCv13**：两版形状的局部法向、凹陷和经度配准并不完全相同；4 μm 位于 Wien 侧，对最热面元和阴影边界非常敏感。
4. **拟合网格抽样**：当前扫描只有 96 个宏面元代表 12,288 个面元；面积守恒不保证温度极值、遮挡和定向热辐射严格守恒。
5. **OVIRS SRF 近似**：公开资料支持高斯线形，但没有取得作者使用的逐像元实验室响应表。本项目以中心 4.00038 μm、FWHM 0.00641 μm 卷积。
6. **全局—粗糙度耦合仍近似**：坑内过程完整迭代；全局发射源目前使用宏观光滑面温度，没有对每个坑微面元和所有远处宏面元做完全耦合迭代。
7. **误差模型不完备**：1° 箱误差包含统计散布、传播误差和 0.5% 下限，却没有作者的系统协方差。约化 χ²=154.24 表明残差由系统误差主导，Δχ²=1 不能作为可信物理误差条。

## 7. GitHub 与 VS Code 上传说明

大型 PDS 数据不进入 Git。`.gitignore` 排除了 19 GB OVIRS FITS、SPICE 二进制核、OTES、视因子缓存、临时文件和大输出；GitHub 中保留下载器，使其他人可以从官方源重建数据。

VS Code 正确操作：

1. 打开 `D:\ChatGPT\热物理反演`；
2. Source Control 中确认不会出现 `data/ovirs/full`；
3. 只暂存源代码、Markdown、小型形状和最终摘要；
4. Commit；
5. 点击 Sync Changes，或在终端执行 `git push`。

如果 VS Code 一直转圈，先执行 `Developer: Reload Window`。如果输出 `RPC failed` 或推送数百 MB，说明历史里已有大二进制文件；`.gitignore` 只影响新文件，不能删除旧提交中的对象。本项目当前历史最大单文件约 57 MB，没有超过 GitHub 的 100 MB 单文件硬限制，认证预检也已成功。

