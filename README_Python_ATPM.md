# Bennu ATPM-like Python 复现

本项目现在包含一个仅依赖 `numpy` 的 Python 热物理正演与反演流程。它在现有单面元脚本基础上增加了三角形形状模型、盘积分、周期稳态一维导热、半球坑微面元粗糙度、4/14 μm 普朗克辐射、χ² 网格和 100 次 bootstrap。

## 直接运行

正式运行（论文 Γ 网格 0–600，步长 10；约需数分钟）：

```powershell
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' .\reproduce_bennu.py
```

快速检查：

```powershell
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' .\reproduce_bennu.py --quick --bootstrap 10
```

测试：

```powershell
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -m unittest -v .\test_bennu_atpm.py
```

若获得任务使用的 Bennu OBJ 网格，可替换内置代理形状：

```powershell
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' .\reproduce_bennu.py --shape 'D:\data\bennu_v13.obj'
```

## 建模过程

1. `make_bennu_proxy_mesh()` 建立带赤道脊和小尺度纵向变化的 224 面元“陀螺形”网格；`load_obj()` 可读取外部 OBJ。
2. 对每一转相计算面元法向与太阳、观测者方向的点积，得到日照和可见投影。
3. 每个面元求解半无限介质的一维导热。频域中的地下导热通量为 `Γ sqrt(i n ω) T_n`，表面满足 `Q_abs = εσT⁴ + q_cond`。频域求解与 Julia 包的时间步进离散求的是同一个偏微分方程和辐射边界，但能直接得到周期稳态。
4. 一个宏观面元分成光滑部分和完全半球坑部分。生产配置使用 8×16=128 个等立体角坑壁面元，并计算坑口射线遮挡、太阳多次散射、视因子热自加热和方向性出射。相对 10×20=200 面元参考，代表性 SPCv14+SPICE 测试的 4 μm 曲线 RMS 差为 0.77%，最大差为 1.55%。
5. 分别用普朗克定律计算 4 μm（OVIRS 约束 Γ）和 14 μm（OTES 相对光变约束粗糙度）的盘平均光谱辐亮度。
6. 对 Γ=0–600、坑覆盖率=0–1 计算 χ²；用固定随机种子的 100 组合成数据重新拟合，给出 bootstrap 标准差。

## 本次运行结果

|量|论文|Python 正式运行|
|---|---:|---:|
|热惯量 Γ|350 ± 20 J m⁻² K⁻¹ s⁻½|350；bootstrap σ=22.16|
|半球坑覆盖率|0.77 ± 0.04|0.77；bootstrap σ=0.0457|
|RMS 坡度|43 ± 1°|43.0°；bootstrap σ=1.29°|
|名义昼夜趋肤深度（`cp=800, rho=2500`）|论文参数导出|1.228 cm|

输出目录 `output/python_bennu/`：

- `summary.json`：参数、结果、物理检查和未包含的过程；
- `best_fit_curves.csv`：4 μm 与 14 μm 的代理观测及最优曲线；
- `chi2_grid.csv`：全部 6161 个 Γ–粗糙度组合；
- `bootstrap.csv`：100 次反演结果；
- `lightcurve_comparison.svg`：曲线对照图。

## 结果边界（很重要）

这里恢复论文数值，证明的是代码链条能够把“形状 → 温度 → 辐亮度 → χ² → 参数和误差”完整跑通。代理观测由论文最优参数生成，噪声幅度按论文公开的不确定度量级校准，因此它不是对论文结论的独立重新测量。

严格复现仍需要论文没有随附件发布的输入：作者的 L3a 逐点辐亮度/误差和欠填充视场经验修正、精确 SPCv13 网格及像元级实验室光谱响应。当前公开数据流程已经接入全部帧的 SPICE SPK/CK/IK 几何、OVIRS 实际视轴与 4 mrad 视场、公开孔径响应近似、SPCv14 全局形状光线追踪/相互辐射，以及半球坑的遮挡、热自加热、多次太阳散射和 thermal beaming。OVIRS 4 μm 响应采用由公开波长采样和仪器论文给出的高斯线形构成的卷积，不冒充未公开的逐像元 SRF。

## 真实 PDS 数据反演

`pds_observations.py` 按 PDS4 XML 定义读取公开数据，而不是从论文图中取点：

- OTES：读取两个 24 MB L2 二进制文件的 2×8,675 条记录；通过 SCLK 与 17,350 条 geometry 逐条匹配；要求 `look_type=data-look`、`bore_flag=1`，并按 OTES SIS 的 quality 位筛选；每个观测日保留 7,839 条记录，在 14 μm 附近平均并分成 96 个自转相位箱。
- OVIRS：读取 L2 FITS 的 20×512 辐亮度、波长、质量和不确定度数组；按 OVIRS SIS 排除空超像元和 outlier；用 2.05–2.15 μm 拟合太阳光谱形状并扣除反射分量。论文拟合使用 3.5–4.0 μm 光谱，而补充图 2b 是其中的 4 μm 通道；因此单色图和当前单色模型使用 3.98–4.02 μm 样本，不能把整个 3.5–4.0 μm 平均值与 4.0 μm 模型比较。项目中现有 2018-11-02/03 的 34,254 个官方校准 FITS（19,138,394,880 字节），全部通过固定大小和 FITS 文件头校验。

正式真实数据反演：

```powershell
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  .\invert_real_bennu.py `
  --gamma-step 10 `
  --fraction-step 0.01 `
  --phase-bins 96 `
  --output .\output\real_bennu
```

快速检查：

```powershell
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  .\invert_real_bennu.py `
  --gamma-step 50 `
  --fraction-step 0.05 `
  --phase-bins 64 `
  --output .\output\real_bennu_quick
```

真实数据输出包括 `real_summary.json`、`real_lightcurves.csv`、`real_chi2_grid.csv` 和 `real_lightcurve_fit.svg`。本次正式网格的最优值落在允许范围的上边界（Γ=600、坑覆盖率=1），约化 χ² 约为 39.45。这不是新的 Bennu 物性测量；它说明真实 OTES 曲线不能由当前 224 面元代理形状、平均观测几何和简化粗糙度模型充分解释。只有补齐完整 OVIRS、自洽的逐面元 SPICE 几何、论文 SPC v13 网格和作者 ATPM 过程后，才可以将真实数据反演值与论文的 Γ=350±20 直接比较。

### 2018-11-02/03 OVIRS 与补充图 2b

官方 PDS 数据入口为 https://sbnarchive.psi.edu/pds4/orex/orex.ovirs/data_calibrated/approach/ 。下载器限定两天的校准 L2 FITS，支持断点续传并在结束时写入下载状态：

~~~powershell
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -u .\download_ovirs_20181102_03.py --workers 24 --output .\data\ovirs\full
~~~

先下载官方最小 SPICE 核集合，再按论文处理并反演（Windows 下从 ASCII 联接 `D:\ATPM` 运行）：

~~~powershell
cd D:\ATPM
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' .\download_spice_kernels.py
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -u .\reproduce_supplementary_fig2b.py --input .\data\ovirs\full --shape .\data\shape\g_12560mm_spc_obj_0000n00000_v014.obj --output .\output\supplementary_fig2b --gamma-min 250 --gamma-max 400 --gamma-step 2 --model-phases 64 --fit-facet-stride 16
~~~

公开 L2 探测器辐亮度在 2018-11-02 的均值为 1.1235×10⁻⁴ W cm⁻² μm⁻¹ sr⁻¹；除以 PDS `FILL_FAC` 后，盘面等效均值为 3.0716×10⁻⁴。二者定义不同，不能交叉比较。现已下载官方 SPCv14 12.56 m 替代形状（6,534 顶点、12,288 面元，PDS km 坐标自动换算为 m），并在坑尺度和全局形状尺度加入相互辐射与射线可见性。全部 34,254 帧通过重建 CK 的姿态、IK 的 OVIRS 视轴/FOV、SPK 轨道与 Bennu PCK 自转逐帧求几何；首帧 CK 视轴与 FITS 头的 RA/Dec 一致到约 4×10⁻⁶ 度。

当前 128 坑面元、96 个均匀抽样全局面元、48 相位的精细扫描得到 Γ=248 J m⁻² K⁻¹ s⁻½，但约化 χ²=154.24，明显不是可接受的统计拟合，也不能替代论文的 350±20。这个偏差主要说明公开 L2→论文 L3a 的欠填充 FOV 经验修正未公开，且拟合仍对全局形状作了 1/128 抽样。详见 `output/supplementary_fig2b/summary.json` 和 `PHYSICS_UPGRADE_REPORT.md`。

SPICE 的概念、核文件职责、公式和字段说明见 `SPICE_GEOMETRY.md`；坑分辨率数据见 `output/crater_resolution_convergence.csv`。

低内存运行可使用：

~~~powershell
& 'C:\Users\松廷\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -u .\reproduce_supplementary_fig2b.py --shape .\data\shape\g_12560mm_spc_obj_0000n00000_v014.obj --gamma-min 0 --gamma-max 1200 --gamma-step 20 --fit-facet-stride 16 --facet-chunk-size 32 --self-heating-iterations 4
~~~

## 与下载的 Julia 包关系

两者一致的部分是面元日照、`k=Γ²/(ρcp)` 对应的一维导热物理、非线性辐射边界、普朗克辐射和盘积分思想。Python 版用频域周期解代替 Julia 的显式/隐式/Crank–Nicolson 时间推进，并实现了本任务需要的分数半球坑、Γ–粗糙度网格及 bootstrap。下载的 Julia 包当前也不等于论文作者的 ATPM：其路线图仍将 roughness-aware problem 列为未实现，而且没有论文的 OVIRS/OTES 似然和任务输入。
