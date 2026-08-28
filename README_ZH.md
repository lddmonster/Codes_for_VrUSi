# Multi-frame Reconstruction via Clutter Filtering for Computational Contrast-Enhanced Ultrasound Microflow Imaging

> 主函数：`PER_L12_64eles_veraRF_recon.m`（本文件夹根目录）
> 本文档包含：文件分类结构、算法解析（流程 / 目标函数 / 维度 / 标准差 / 缩放 / τ / TwIST 细节）、使用方法
> 符号约定：λ = c/f₀ ≈ 197.1 μm（用于换算实际尺寸）；T = 帧数 = 1000

---

## 目录

1. [文件结构（分类存放）](#1-文件结构分类存放)
2. [整体流程（伪代码）](#2-整体流程伪代码)
3. [精确目标函数](#3-精确目标函数)
4. [参数维度](#4-参数维度)
5. [正则项中的标准差计算](#5-正则项中的标准差计算)
6. [每个正则项的缩放方式](#6-每个正则项的缩放方式)
7. [惩罚参数 τ 的使用](#7-惩罚参数-τ-的使用)
8. [TwIST 实现细节](#8-twist-实现细节)
9. [使用方法](#9-使用方法)

---

## 1. 文件结构（分类存放）

```
PER_L12_64eles_veraRF_recon_deps\
├── README.md / README_ZH.md     本文档（英文 / 中文）
├── PER_L12_64eles_veraRF_recon.m               主函数（入口）
├── 01_main_pipeline\            主流程脚本（按执行顺序）
│   ├── L12_64eles_veraRF_recon_continue1.m    阵元/波形/成像网格参数定义
│   ├── L12_64eles_veraRF_recon_continue2.m    RF 数据装载 → 观测矩阵 Y
│   ├── L12_64eles_veraRF_recon_continue3.m    脉冲响应装载 → 系统矩阵 R
│   └── TWIST_vera_L12_64eles.m                组装 TwIST 优化问题并求解
├── 02_optimization\             优化求解器与算子
│   ├── TwIST.m                                 两步迭代收缩阈值求解器
│   ├── soft2.m                                 复数软阈值（临近算子 Psi）
│   └── weighted_L1.m                           正则项 Φ（std-L1 + 逐帧平均 L1）
├── 03_output_io\                结果输出
│   ├── tifwrite.m                              3D tif 写入封装
│   └── Fast_Tiff_Write.m                       tif 写入类
├── 04_Field_II\                 Field II 仿真库（本流程仅用初始化/采样率）
│   ├── field_init.m
│   ├── set_sampling.m
│   └── Mat_field.mexw64                        Field II 的 C 接口
└── 05_MUST\                     MUST 工具箱
    └── getparam.m                              探头参数 getparam('L12-3v')
```

| 文件 | 原路径（相对 `K:\scui code copy ver-20250602\fieldii\pulse-echo impulse response\`） |
|---|---|
| 01_main_pipeline\*.m | `custom10M\`（weighted_L1、TWIST_vera 为本流程修改版） |
| 02_optimization\TwIST.m | `TwIST_v2\TwIST_v2\`（自包含） |
| 02_optimization\soft2.m | `custom10M\` |
| 02_optimization\weighted_L1.m | `.\`（`MPI\` 下有一份旧版，Nmpz=325、返回 y1，与本流程 Nmpz=285 不符，勿混用） |
| 03_output_io\tifwrite.m | `.\`（根目录 `K:\scui code copy ver-20250602\` 下另有支持 4D 的版本，写 3D 数据两者等效） |
| 03_output_io\Fast_Tiff_Write.m | `.\`（各处副本内容相同） |
| 04_Field_II\* | `Field_II_windows\`（完整 Field II 库仍在原处） |
| 05_MUST\getparam.m | `MUST\` |

---

## 2. 整体流程（伪代码）

```
输入: ① Verasonics 采集的 RF 数据集（1000 帧、发射角#3、64 通道）
      ② 预先算好的 Field II 脉冲响应库（64 阵元 × 36480 像素）

── 阶段1 系统定义 (continue1) ──
c=1540, f0=7.813MHz, fs=31.25MHz, 64阵元, pitch=0.2mm
探头参数 ← getparam('L12-3v')                       # MUST 工具箱
加载实测 impulse response (tv.mat, 降采样+归一化)      # 接收/发射特性
加载实测发射波形 (TW_L12_067DC_2cycles.mat)
定义成像网格: 128(横向) × 285(深度) = 36480 像素
             横向 ±32.73λ (≈±6.45mm)，深度 10.15λ~152.2λ (≈2~30mm)

── 阶段2 观测矩阵 (continue2) ──
for 帧 f = 1..1000:
    取角#3、通道65:128 的 RF 段 (len×64), NaN→0
    映射到以 f0 为中心的 1200 点子带 → 展开成 76800×1 列向量
堆叠: Y = rf_all ∈ ℂ^(76800×1000)        # 76800 = 1200×64

── 阶段3 系统矩阵 (continue3) ──
for 阵元 e = 1..64:
    载入该阵元对全部 36480 像素的脉冲响应 (36480 列)
    同一子带截取 → 1200×36480
垂直堆叠 64 阵元: R ∈ ℂ^(76800×36480)    # 与 Y 同一观测空间

── 阶段4 TwIST 稀疏求解 (TWIST_vera_L12_64eles) ──
τ  ← 0.02 · max|R'Y|
X₀ ← 0 (36480×1000)
X* ← TwIST(Y, R, τ;  Psi=soft2, Phi=weighted_L1, λ₁=1e-4,
           Monotone=1, StopCriterion=1, tol=1e-2, Maxiter=50)

── 阶段5 输出 (主函数) ──
for f = 1..1000: img_f = reshape(|X*(:,f)|, 285, 128), 线性缩放到 [0,255]
写成 3D tif (285×128×1000), 并保存 X*.mat
```

**调用链**：主函数 → continue1（参数）/ continue2（RF→Y）/ continue3（脉冲响应→R）→ TWIST_vera（组装优化问题）→ TwIST.m（求解器，内部用 soft2 作临近算子、weighted_L1 作目标评估）→ 主函数后处理输出。

---

## 3. 精确目标函数

`TwIST.m` 中每次迭代评估的目标为 `f = 0.5*resid'resid + tau*phi(x)`，本代码实际求解的问题是：

$$\min_{X\in\mathbb{C}^{36480\times 1000}}\; F(X)=\tfrac{1}{2}\big\|Y-RX\big\|_F^2+\tau\,\Phi(X)$$

$$\Phi(X)=\underbrace{\big\|\sigma(X)\big\|_1}_{\text{时间标准差图稀疏}}\;+\;\underbrace{\tfrac{1}{T}\big\|X(:)\big\|_1}_{\text{逐帧平均的全局 L1}},\qquad \tau = 0.02\times\max\big|R'Y\big|$$

其中逐像素时间标准差（T = 1000 帧，按 N−1 无偏归一化，复数按模长偏差计算）：

$$\sigma_p=\sqrt{\tfrac{1}{T-1}\sum_{t=1}^{T}\big|X_{p,t}-\bar x_p\big|^2},\qquad \bar x_p=\tfrac{1}{T}\sum_{t=1}^{T}X_{p,t}$$

展开即：

$$\Phi(X)=\sum_{p=1}^{36480}\sigma_p\;+\;\frac{1}{T}\sum_{p=1}^{36480}\sum_{t=1}^{T}\big|X_{p,t}\big|$$

第一项约束**时间波动（血管结构）稀疏**，第二项约束**时空幅度（微血流轨迹）稀疏**；两项均为逐像素/逐元素量级，直接相加。

（`'Lambda'` 传入的 λ₁=1e-4 不进入目标函数，仅用于计算 TwIST 迭代步长参数 α、β，见第 8 节。）

---

## 4. 参数维度

| 符号/变量 | 维度 | 说明 |
|---|---|---|
| `Y` = `rf_all(:,1:1000)` = y | **76800 × 1000**，复数 | 观测/RF 矩阵；76800 = 1200 子带点 × 64 接收阵元；列 = 慢时间帧 |
| `R` = `simData_all` | **76800 × 36480**，复数 | 系统矩阵；行同 Y 的观测空间，列 = 像素 |
| `x_twist` = X* | **36480 × 1000**，复数 | 优化变量/重建图像序列；36480 = 285 深度 × 128 横向（像素按 x 块、每块 285 个深度排序） |
| `initarray` | 36480 × 1000 | 全零初始化 |
| `hRt(y)` = R'Y, `grad` | 36480 × 1000 | 梯度/反投影 |
| `resid` = Y−RX | 76800 × 1000 | 残差 |
| `tau` | 标量 0.02 | 惩罚系数因子，实际惩罚参数 = `tau·max(max(abs(R'Y)))`（代码也支持 size(x) 的逐像素 τ，此处未用） |
| `lambda1` | 标量 1e-4 | 只影响 α、β，不在目标中 |
| `bm_kk` | 285 × 128 | 单帧幅度图 |
| `temp_bm_kk` | 285 × 128 × 1000, single | 输出 tif 体数据 |
| weighted_L1 内 `xstack` | 285 × 128 × 1000 | `reshape(x,285,128,[])`，第 3 维 = 帧 |
| weighted_L1 内 `x_mip2` = σ | 285 × 128 | 帧方向 std 后的时间波动图 |
| `y2` = ‖σ‖₁, `y3` = ‖X(:)‖₁, `y3_scale` = y3/T | 标量 | 正则项中间量 |

物理网格：横向 ±6.45 mm、深约 2–30 mm；每帧 `RcvDataFrame` 取通道 `65:128`（64 通道）、角度 `3:3`（只用 1 个发射角）。

---

## 5. 正则项中的标准差计算

`02_optimization\weighted_L1.m`：

```matlab
xstack  = reshape(x, 285, 128, []);   % 像素 × 横向 × 帧(慢时间)
x_mip2  = std(xstack, 0, 3);          % 沿第 3 维(帧/慢时间轴)逐像素标准差
y2      = norm(x_mip2, 1);            % ‖σ‖₁
```

- `std(xstack,0,3)`：**沿第 3 维（帧/慢时间轴）逐像素计算时间标准差**，权重标志 0 表示按 **N−1（=999）归一化**的无偏估计：

$$\sigma_p=\sqrt{\tfrac{1}{T-1}\sum_{t=1}^{T}\big|x_{p,t}-\bar x_p\big|^2}$$

- 复数数据下 MATLAB 的 std 以**偏差的模长**参与求和，并对复数均值 x̄_p 取模长偏差，因此 σ_p 是复像素时间序列的幅度波动量。
- 得到的 285×128 的 σ 图本质上是**能量多普勒式的像素时间波动图**（对血流/血管敏感）；`y2 = ‖σ‖₁` 即"时间波动稀疏"约束，是当前正则 Φ 的**第一项（生效）**。

---

## 6. 每个正则项的缩放方式

`02_optimization\weighted_L1.m`：

| 项 | 表达式 | 缩放方式 | 状态 |
|---|---|---|---|
| **y2** | $\|\sigma(X)\|_1$（时间标准差图 L1，约束血管结构稀疏） | 无缩放 | **生效**（第一项） |
| y3 | $\|X(:)\|_1$（时空全局 L1，约束微血流轨迹稀疏） | 无缩放 | 中间量 |
| **y3_scale** | $y_3/T$，T = `size(x,2)` = 1000 | **除以帧数** | **生效**（第二项） |
| **Φ = y** | $y_2 + y_3/T$ | 两项直接相加 | 最终正则 |

**缩放的作用**：$\|X(:)\|_1$ 随帧数 T 线性增长（重建越多帧该项越大），而 $\|\sigma\|_1$ 大体不随 T 增长（同幅度时间序列的 std 与帧数无关）。将 y3 除以 T 后，两项都折算到"单帧量纲"，组合权重不再随帧数变化——即在"结构稀疏"与"轨迹稀疏"之间取得与 T 无关的平衡。

---

## 7. 惩罚参数 τ 的使用

**取值（`01_main_pipeline\TWIST_vera_L12_64eles.m`）**：

```matlab
tau = 0.02;                                        % 惩罚系数因子（命名常量）
TwIST(y, hR, tau.*max(max(abs(hRt(y)))), ...)      % 实际惩罚参数
```

$$\tau_{\text{eff}} = 0.02\times\max\big|R'Y\big|$$

即匹配滤波/反投影峰值的 **2%**：0.02 作为独立因子抽出，与数据幅度峰值 max|R'Y| 在调用处相乘，效果仍是随数据幅度自适应的无量纲选择（数据整体放大 k 倍，τ_eff 也放大 k 倍）。max|A'y| 正是 LASSO 零解阈值——τ 若 ≥ max|A'y| 解恒为零，取 2% 保持在零解阈值之内。（TwIST 内部的零解快速返回分支只在默认 soft/L1 时启用，本代码传入自定义 Psi，该分支被跳过。）

**在算法中的三个作用点**：

1. 目标函数：`f = 0.5‖resid‖² + tau*Phi(x)`（停机判据 1 与单调性检验都基于它）；
2. 临近算子阈值：IST 步 `x = Psi(xm1 + grad/max_svd, tau/max_svd)`——阈值随自适应步长按 **τ/max_svd** 同步缩放；
3. 不进入梯度：`grad = AT(resid)`，与 τ 无关。

τ 也支持逐像素矩阵形式（`tau` 与 x 同维即可），本代码用标量。

---

## 8. TwIST 实现细节

### 8.1 迭代结构（TwIST.m:481-529）——两步迭代收缩阈值

```
初始化: xm1 = xm2 = X₀ = 0;  max_svd = 1;  Aty = A'y 预计算一次
每轮外层迭代:
    grad ← AT(resid)                      # 残差反投影, 内层循环中不更新
    内层(回溯)循环:
        x ← Psi(xm1 + grad/max_svd, tau/max_svd)        # IST 步
        若已做 ≥2 次 IST 或已进入 TwIST 阶段:
            [sparse=1] mask = (x≠0); xm1,xm2 ← xm1.*mask, xm2.*mask   # 现在为零→过去置零
            xm2 ← (α−β)·xm1 + (1−α)·xm2 + β·x          # 两步外推
            f ← 0.5‖y−A xm2‖² + τΦ(xm2)
            若 f 上升且 Monotone=1: 退回 IST (TwIST_iters=0)
            否则接受 xm2, TwIST_iters++, break
        否则(前两步纯 IST 预热):
            f ← 0.5‖y−A x‖² + τΦ(x)
            若 f 上升: max_svd ← 10·max_svd (缩小步长, 重置计数)   # 处理 ‖A‖²>1 的发散风险
            否则 break
    xm2 ← xm1;  xm1 ← x
    停机: criterion = |f−f_prev|/f_prev (StopCriterion=1)
          继续 = (iter≤50) 且 (criterion>1e-2)；至少跑 5 次迭代
```

### 8.2 步长参数（由 λ₁=1e-4、λ_N=1 计算，TwIST.m:320-326）

- ρ₀ = (1−λ₁/λ_N)/(1+λ₁/λ_N) ≈ 0.99980
- **α = 2/(1+√(1−ρ₀²)) ≈ 1.9608**
- **β = 2α/(λ₁+λ_N) ≈ 3.9212**
- 两步外推系数实际为 (α−β) ≈ −1.96、(1−α) ≈ −0.96、β ≈ 3.92（典型的欠阻尼外推）

### 8.3 其他要点

- **max_svd 自适应**：初始 1（假设 ‖A‖₂²≤1，对未归一化的 R 并不成立）；IST 阶段单调性被破坏时 ×10 回退；每 10000 次 TwIST 迭代 ×0.9（缓慢放大步长）。IST 临近步的阈值 τ/max_svd 与步长 1/max_svd 同步。
- **soft2 即复数软阈值**：`y = max(|x|−T,0); y = y./(y+T).*x`，代数上恒等于 `x·(1−T/|x|)₊`，规避了除零，对复数按模长收缩、保持相位。
- **Phi 与 Psi 分工**：`weighted_L1`（Phi）只用于目标值评估（停机+单调性检验）；真正塑造解的是 `soft2`（Psi）。注意 soft2 是 τ‖·‖₁ 的精确临近算子，与含帧间 std 耦合的 Φ 并不严格匹配——迭代中的收缩行为等价于 L1 软阈值，Φ（含 std 项）通过目标值影响停机与单调性回退，从而影响全部 1000 帧的迭代进程（梯度与矩阵乘仍是逐列作用）。
- **Debias=0**：不做 CG 去偏阶段；输出即主算法解。
- **Verbose=1**：每迭代打印 objective/非零元数。
- 停机判据是**全部 1000 帧的总目标**的相对变化 <1%，个别小幅度帧的收敛由整体主导。

---

## 9. 使用方法

把本文件夹（含全部子文件夹）加入 MATLAB path，并保证数据文件路径可解析：

```matlab
addpath(genpath('...\PER_L12_64eles_veraRF_recon_deps'));            % genpath 自动包含 5 个子文件夹
addpath(genpath('K:\scui code copy ver-20250602\fieldii\pulse-echo impulse response\MPI'));  % tv.mat / TW_L12_067DC_2cycles.mat
run PER_L12_64eles_veraRF_recon.m
```

运行还需的外部数据（非函数，未随文件夹分发）：

| 数据 | 用途 | 位置 |
|---|---|---|
| `tv.mat` | 实测 impulse response（continue1） | `..\L12hydrophone\`、`..\MPI\`（靠 MATLAB path 解析） |
| `TW_L12_067DC_2cycles.mat` | 实测发射波形（continue1） | `..\MPI\` |
| Verasonics RF 数据集 | 观测数据（continue2） | `O:\OPTI\...`、`I:\OPTI\...`（绝对路径） |
| 仿真脉冲响应矩阵 | 系统矩阵（continue3） | `I:\OPTI\Normc_normalized_...`（绝对路径） |

MATLAB 内置/工具箱函数（`rescale`、`downsample` 等）无需复制。
