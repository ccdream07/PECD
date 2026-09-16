# PECD 数据处理系统方案

## 1. 目标与适用范围

本方案面向 Andor 相机采集的 VMI 光电子图像，以及经过 pBasex/Basex Abel 反演的结果。目标是建立一条可重复、可追溯的处理链：

```text
Andor 原始图像
  -> 暗场/背景/坏点处理
  -> LCP/RCP 配对与强度归一化
  -> VMI 中心和方向校正
  -> Abel 反演（pBasex 或 pyAbel）
  -> 半径到电子动能标定
  -> PES、PAD 和 Legendre 系数
  -> PECD、误差和质量控制
  -> conformer/对映体分析与报告
```

当前 `test` 目录中的文件是一个样例，不应把文件名 `200-130-60-2` 当作完整实验元数据。

## 2. 样例文件的初步判断

| 文件 | 当前尺寸 | 用途判断 |
|---|---:|---|
| `200-130-60-2.dat` | `1024 x 1024` | 原始二维相机图像 |
| `200-130-60-2_reconst.dat` | `1011 x 1015` | Abel/Basex 重构图 |
| `200-130-60-2_reconst.bmp` | `1091 x 1011` | 重构结果的可视化图 |
| `200-130-60-2_polar.dat` | `159840 x 3` | 半径、角度、强度；约为 `999 x 160` 个点 |
| `200-130-60-2_speed.dat` | `1000 x 2` | 半径和径向积分强度 |
| `200-130-60-2_aniso.dat` | `999 x 4` | 半径及三个各向异性列；列含义需由 Basex 配置确认 |

样例中半径大约为 `0.5--499.5`，角度约为 `9.5--89` 度。`aniso.dat` 在低信号大半径处出现很大的系数，说明正式程序必须按有效计数或信噪比过滤，不能直接使用所有半径点。

## 3. 文献中的核心公式

### 3.1 光电子角分布

对于圆偏振光，光电子角分布写为：

\[
I^{\{p\}}(\theta)=
\frac{I_{\mathrm{tot}}^{\{p\}}}{4\pi}
\left[1+\sum_{i=1}^{2n}\beta_i^{\{p\}}P_i(\cos\theta)\right],
\]

其中 (p=+1,-1) 表示两种相反圆偏振，(P_i) 是 Legendre 多项式，\(\beta_i\) 是角分布系数，\(n\) 是吸收的光子数。

二光子电离时通常拟合到 (P_4)：

\[
I(\theta)=I_0\left[1+\beta_1P_1+\beta_2P_2+\beta_3P_3+\beta_4P_4\right].
\]

奇数项 \(\beta_1,\beta_3\) 负责前后不对称，偶数项描述普通 PAD 形状。对于非手性分子或线偏振光，奇数项应为零；反转圆偏振手性时，奇数项应反号。

### 3.2 二光子 PECD

Rouquet 等文献使用：

\[
\boxed{\mathrm{PECD}=2\beta_1-\frac{1}{2}\beta_3}
\]

这里的 \(\beta_1,\beta_3\) 应明确是左右圆偏振差分约定下的“奇数项幅度”。推荐在程序中定义：

\[
\beta_l^{\mathrm{odd}}=
\frac{\beta_{l,\mathrm{LCP}}-\beta_{l,\mathrm{RCP}}}{2},
\]

再使用：

\[
\mathrm{PECD}=2\beta_1^{\mathrm{odd}}-\frac12\beta_3^{\mathrm{odd}}.
\]

如果程序直接拟合 LCP-RCP 差图得到 \(\Delta\beta_l=\beta_{l,\mathrm{LCP}}-\beta_{l,\mathrm{RCP}}\)，等价公式是：

\[
\boxed{\mathrm{PECD}=\Delta\beta_1-\frac14\Delta\beta_3}.
\]

两种写法只能选一种，避免重复乘以 2。

### 3.3 半球积分 PECD

对每种 helicity 计算前向和后向积分：

\[
A_h=\frac{Y_{\mathrm{fwd}}^h-Y_{\mathrm{bwd}}^h}
{Y_{\mathrm{fwd}}^h+Y_{\mathrm{bwd}}^h}.
\]

文献给出的半球法为：

\[
\boxed{\mathrm{PECD}=2(A_{\mathrm{LCP}}-A_{\mathrm{RCP}})}.
\]

该结果应作为 PAD/Abel 结果的独立交叉检查。半球积分必须使用相同的半径或电子动能 ROI，并明确光传播方向和正负号。

### 3.4 半径到电子动能

VMI 中常用标定形式为：

\[
E_{\mathrm{kin}}=k(r-r_0)^2,
\]

必要时使用：

\[
E_{\mathrm{kin}}=a_0+a_1(r-r_0)^2+a_2(r-r_0)^4.
\]

参数必须由已知能级、已知光电子峰或独立标定实验给出，不能从当前样例文件名推断。

## 4. 推荐处理步骤

### 4.1 数据登记

每次采集必须保存：

- 原始 Andor 文件和文件格式；
- LCP/RCP 标签及采集顺序；
- 样品、对映体、构象和共振波长；
- 激发光/电离光波长、光子数、脉冲能量和偏振；
- 相机曝光、增益、位深和帧数；
- VMI 电极电压、图像中心和传播方向；
- 暗场、背景和饱和掩膜；
- 半径到能量的标定参数；
- pBasex/pyAbel 方法、基函数阶数、正则化和版本。

建议每个实验目录有一个 `metadata.yaml`，原始数据只读保存，所有处理结果写入新的处理目录。

### 4.2 原始图像预处理

\[
I_{\mathrm{clean}}=I_{\mathrm{raw}}-I_{\mathrm{dark}}-I_{\mathrm{background}}.
\]

同时执行坏点/热像素处理、饱和检查、负值策略、探测器响应校正和激光脉冲能量归一化。所有阈值和掩膜都要写入处理记录。

### 4.3 LCP/RCP 配对

对交替采集的左右 helicity：

\[
I_{\mathrm{sum}}=I_{\mathrm{LCP}}+I_{\mathrm{RCP}},\qquad
I_{\mathrm{diff}}=I_{\mathrm{LCP}}-I_{\mathrm{RCP}}.
\]

`sum` 用于 PES/PAD，`diff` 用于 PECD。应检查长期漂移、两种 helicity 的总计数和配对顺序。

### 4.4 中心、方向和 ROI

\[
r=\sqrt{(x-x_0)^2+(y-y_0)^2},\qquad
\theta=\arccos\left(\frac{y-y_0}{r}\right).
\]

必须校准：图像中心、光传播方向、前向/后向定义、中心低能掩膜、最大有效半径和饱和区域。PECD 正负号应由已知样品或已知 helicity 实验确认。

### 4.5 Abel 反演与 PAD

对 LCP、RCP 或其 sum/diff 进行 pBasex/pyAbel 反演，保存重构残差。对每个半径/能量 bin 拟合：

\[
I_h(r,\theta)=I_{0,h}(r)\left[1+\sum_l\beta_{l,h}(r)P_l(\cos\theta)\right].
\]

输出 \(I_0\)、\(\beta_1\)--\(\beta_4\)、PECD、拟合残差、协方差、有效计数和信噪比。

### 4.6 误差

计数近似为 Poisson 分布时，对

\[
A=\frac{F-B}{F+B}
\]

可使用：

\[
\mathrm{Var}(A)\approx
\frac{4B^2F+4F^2B}{(F+B)^4}.
\]

PECD 误差近似为：

\[
\mathrm{Var}(\mathrm{PECD})=
4[\mathrm{Var}(A_{\mathrm{LCP}})+\mathrm{Var}(A_{\mathrm{RCP}})].
\]

正式结果还应加入配对 bootstrap、激光漂移、中心位置、ROI、能量标定和拟合协方差造成的不确定度。

## 5. Conformer 和对映体分析

窄带 (S_1\leftarrow S_0) 共振可选择特定 conformer 或振动态。每个选择条件应独立输出：

\[
\mathrm{PECD}_c(E_{\mathrm{kin}}).
\]

多个构象混合时，不能对 PECD 做无权平均。应按信号强度加权：

\[
I_{\mathrm{tot}}(E)=\sum_c w_cI_c(E),
\]

\[
\boxed{
\mathrm{PECD}_{\mathrm{tot}}(E)=
\frac{\sum_c w_cI_c(E)\mathrm{PECD}_c(E)}
{\sum_c w_cI_c(E)}}.
\]

单一 conformer 且响应线性时，可用：

\[
\mathrm{PECD}_{\mathrm{obs}}=ee\cdot\mathrm{PECD}_{\mathrm{pure}},
\qquad
ee=\frac{\mathrm{PECD}_{\mathrm{obs}}}{\mathrm{PECD}_{\mathrm{pure}}}.
\]

多构象、不同电子能量或不同跃迁重叠时，应使用显式混合模型。

## 6. 建议编写的程序

1. `data_schema.py`：实验元数据和版本定义。
2. `andor_io.py`：读取 `.dat`、`.sif`、`.tif` 等 Andor 数据。
3. `preprocess.py`：暗场、背景、坏点、饱和和归一化。
4. `image_registration.py`：中心、方向、镜像和坐标校正。
5. `helicity_pairing.py`：LCP/RCP 配对和漂移检查。
6. `vmi_calibration.py`：半径到电子动能标定。
7. `abel_inversion.py`：pBasex/pyAbel 反演。
8. `pbasex_reader.py`：读取和验证 `speed/polar/aniso/reconst`。
9. `pad_fit.py`：Legendre 系数拟合。
10. `pecd_analysis.py`：半球法和 PAD 法 PECD。
11. `uncertainty.py`：Poisson、bootstrap 和系统误差。
12. `conformer_analysis.py`：构象分组和混合模型。
13. `quality_control.py`：反号、racemate、残差、漂移和低计数检查。
14. `report.py`：图表、公式、参数和结果报告。
15. `cli.py`：`ingest`、`preprocess`、`invert`、`calibrate`、`fit-pad`、`calculate-pecd`、`report` 命令。

## 7. 第一版最小可用系统

第一阶段建议实现：读取当前样例、生成 sum/diff、指定或读取图像中心、半球积分、读取 Basex 输出、半径到能量转换、能量分辨 PECD、Legendre 系数曲线和 QC 报告。自动 Abel 反演、构象混合拟合和 ee 分析放在第二阶段。

## 8. 质量控制标准

- racemic 样品的 PECD 应接近零；
- 反转 helicity 后 PECD 应反号；
- 反转对映体后 PECD 应近似反号；
- LCP/RCP 的总强度不能存在未解释的大漂移；
- 低计数半径不得产生被误认为物理信号的巨大 \(\beta_l\)；
- Abel 重构残差、拟合残差和有效计数必须随结果保存；
- 半球法与 PAD 法应在相同 ROI 下相互验证；
- `aniso.dat` 的列定义必须由 Basex 配置或输出说明确认。

## 9. 当前必须补充的实验信息

在正式开发前，需要确定：LCP/RCP 文件如何命名和标记、原始文件是否为累计图还是逐帧图、`aniso.dat` 三列的确切含义、图像中心、VMI 能量标定、光传播方向、光子数、偏振纯度、背景采集方式以及 pBasex 参数文件。

