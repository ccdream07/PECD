# 苯甲醚参考样品与 PyAbel 基线分析

苯甲醚是非手性分子。在理想圆偏振实验中，它不应出现真实 PECD。它适合用于检查图像中心、背景扣除、左右圆偏振配对和分析流程是否产生假手性信号。

这不代表它的光电子角分布必定各向同性。\(\beta_2\) 和 \(\beta_4\) 是偶数阶 PAD 系数，能够非零；径向能谱也取决于分子的电离和实验条件。

外消旋手性样品的理想总信号同样没有奇数阶手性项，但它的 PES、\(\beta_2\)、\(\beta_4\) 不必与苯甲醚相同。可比较的是：

\[
\mathrm{PECD}_{\mathrm{anisole}}(E) \approx 0,
\qquad
\mathrm{PECD}_{\mathrm{racemate}}(E) \approx 0.
\]

两者都不是要求 \(\beta_2=\beta_4=0\) 的“各向同性”条件。

## 单图像基线程序

`scripts/analyze_vmi_sample.py` 读取单个二维 `.dat` 图像，使用 PyAbel 的 rBasex：

- 自动估计图像中心，或接受手动中心；
- 生成 rBasex Abel 重构图；
- 导出径向强度与 \(\beta_1\)--\(\beta_6\)；
- 生成径向谱、系数曲线和现有 Basex `aniso.dat` 的诊断图；
- 对低于最大径向强度 1% 的区域只保留曲线，不将其纳入 QC 汇总。

运行当前样例：

```powershell
python scripts\analyze_vmi_sample.py test\200-130-60-2.dat
```

本样例自动找到的中心为 `(row=505, col=507)`。需要人工固定中心时：

```powershell
python scripts\analyze_vmi_sample.py test\200-130-60-2.dat --center-row 505 --center-col 507
```

输出写入 `test/200-130-60-2_analysis/`，包括：

- `raw_image.png`
- `rbasex_reconstruction.png`
- `radial_spectrum.png`
- `anisotropy_qc.png`
- `vendor_aniso.png`
- `rbasex_radial_anisotropy.csv`
- `single_image_odd_component.csv`
- `summary.json`

## Basex 对照验证（2026-09-17）

之前的 `vendor_aniso.png` 把旧 Basex 的三列仅标作 column 1--3，并与 PyAbel 的
`beta_1`--`beta_3` 放在同一张图中。这种比较是不成立的，因为两边的列不是同一阶
Legendre 系数。

对 `200-130-60-2_polar.dat` 按每个半径拟合

\[
I(r,\theta)=c_0(r)+c_2(r)P_2(\cos\theta)+c_4(r)P_4(\cos\theta)+c_6(r)P_6(\cos\theta),
\]

并计算 `beta_l=c_l/c_0` 后，结果与 `200-130-60-2_aniso.dat` 的三列逐点吻合：
中位绝对误差约为 `2.4e-6`、`2.5e-6`、`3.1e-6`，相关系数分别大于
`0.9999998`。因此该 Basex 文件的列定义为

```text
radius_px, beta_2, beta_4, beta_6
```

程序现在默认使用 `order=6`，并在 `vendor_aniso.png` 中把 Basex 的
`beta_2/beta_4/beta_6` 与 PyAbel 的同阶系数叠加比较。比较时必须统一中心、半径坐标、
背景处理、阶数和正则化。当前样例在高信号半径内，\(\beta_2\) 的趋势相关性约为 0.987，
而 \(\beta_4\)、\(\beta_6\) 的相关性约为 0.72、0.82；这说明主要结构有一致性，但高阶项
仍受算法参数和噪声影响。大半径和零信号区的巨大尖峰是系数除以很小的 `c_0` 后放大的
反演噪声，应按强度/SNR 阈值剔除。

这验证了旧 Basex 输出的列含义，但没有证明 PyAbel 与旧 Basex 在所有参数下完全等价。
正式 PECD 仍需成对的 LCP/RCP 图像；单张图的奇数项只能作为对称性和数据质量检查。

## 重要限制

单张图像不能计算定量 PECD。程序输出的

\[
2\beta_1-\frac12\beta_3
\]

只用于单图像前后对称性质量检查，不能报告为 PECD。定量 PECD 需要同一实验条件、明确 helicity 标签的 LCP/RCP 图像对。将来应计算：

\[
\beta_l^{\mathrm{odd}}=
\frac{\beta_{l,\mathrm{LCP}}-\beta_{l,\mathrm{RCP}}}{2},
\]

随后：

\[
\mathrm{PECD}=2\beta_1^{\mathrm{odd}}-\frac12\beta_3^{\mathrm{odd}}.
\]
