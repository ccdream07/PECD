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
- 导出径向强度与 \(\beta_1\)--\(\beta_4\)；
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

