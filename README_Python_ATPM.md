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
4. 一个宏观面元分成光滑部分和完全半球坑部分。坑内由 18 个微面元积分；每一转相按坑口投影重新归一化太阳输入和出射投影，保证能量及投影面积守恒。
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

严格复现仍需要论文没有随附件发布的输入：OVIRS/OTES 逐点辐亮度和误差、每个样本的观测几何/SPICE、当时使用的 SPC v13 约 12 m 网格、仪器波段响应，以及作者定制 ATPM 的阴影/自加热设置。当前模型也明确省略全局光线追踪阴影、迭代自加热、多次散射和仪器响应。拿到这些输入后，可保留导热和反演接口，替换 `Mesh`、方向数组和观测似然。

## 与下载的 Julia 包关系

两者一致的部分是面元日照、`k=Γ²/(ρcp)` 对应的一维导热物理、非线性辐射边界、普朗克辐射和盘积分思想。Python 版用频域周期解代替 Julia 的显式/隐式/Crank–Nicolson 时间推进，并实现了本任务需要的分数半球坑、Γ–粗糙度网格及 bootstrap。下载的 Julia 包当前也不等于论文作者的 ATPM：其路线图仍将 roughness-aware problem 列为未实现，而且没有论文的 OVIRS/OTES 似然和任务输入。
