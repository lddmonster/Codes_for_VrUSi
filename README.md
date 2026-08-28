# Multi-frame Reconstruction via Clutter Filtering for Computational Contrast-Enhanced Ultrasound Microflow Imaging

---

## Table of Contents

1. [File Structure (Categorized)](#1-file-structure-categorized)
2. [Overall Pipeline (Pseudocode)](#2-overall-pipeline-pseudocode)
3. [Exact Objective Function](#3-exact-objective-function)
4. [Parameter Dimensions](#4-parameter-dimensions)
5. [Standard-Deviation Computation in the Regularizer](#5-standard-deviation-computation-in-the-regularizer)
6. [Scaling of Each Regularization Term](#6-scaling-of-each-regularization-term)
7. [Use of the Penalty Parameter τ](#7-use-of-the-penalty-parameter-τ)
8. [TwIST Implementation Details](#8-twist-implementation-details)
9. [Usage](#9-usage)

---

## 1. File Structure (Categorized)

```
PER_L12_64eles_veraRF_recon_deps\
├── README.md / README_ZH.md     This document (English / Chinese)
├── PER_L12_64eles_veraRF_recon.m               Main script (entry point)
├── 01_main_pipeline\            Main pipeline scripts (in execution order)
│   ├── L12_64eles_veraRF_recon_continue1.m    Array / waveform / imaging-grid parameter definition
│   ├── L12_64eles_veraRF_recon_continue2.m    RF data loading → observation matrix Y
│   ├── L12_64eles_veraRF_recon_continue3.m    Impulse-response loading → system matrix R
│   └── TWIST_vera_L12_64eles.m                Assembles and solves the TwIST optimization problem
├── 02_optimization\             Optimizer and operators
│   ├── TwIST.m                                 Two-step iterative shrinkage/thresholding solver
│   ├── soft2.m                                 Complex soft-thresholding (proximal operator Psi)
│   └── weighted_L1.m                           Regularizer Φ (std-L1 + per-frame-averaged L1)
├── 03_output_io\                Result output
│   ├── tifwrite.m                              3D tif writing wrapper
│   └── Fast_Tiff_Write.m                       tif writing class
├── 04_Field_II\                 Field II simulation library (only init/sampling used in this pipeline)
│   ├── field_init.m
│   ├── set_sampling.m
│   └── Mat_field.mexw64                        Field II C-interface MEX file
└── 05_MUST\                     MUST toolbox
    └── getparam.m                              Probe parameters, getparam('L12-3v')
```

| File | Original path (relative to `K:\scui code copy ver-20250602\fieldii\pulse-echo impulse response\`) |
|---|---|
| 01_main_pipeline\*.m | `custom10M\` (weighted_L1 and TWIST_vera are the versions modified for this pipeline) |
| 02_optimization\TwIST.m | `TwIST_v2\TwIST_v2\` (self-contained) |
| 02_optimization\soft2.m | `custom10M\` |
| 02_optimization\weighted_L1.m | `.\` (an older version exists in `MPI\` with Nmpz=325 and returning y1 — incompatible with this pipeline's Nmpz=285; do not mix them up) |
| 03_output_io\tifwrite.m | `.\` (a 4D-capable variant exists at `K:\scui code copy ver-20250602\`; both are equivalent for 3D data) |
| 03_output_io\Fast_Tiff_Write.m | `.\` (all copies identical) |
| 04_Field_II\* | `Field_II_windows\` (the complete Field II library remains there) |
| 05_MUST\getparam.m | `MUST\` |

---

## 2. Overall Pipeline (Pseudocode)

```
Inputs: ① Verasonics-acquired RF dataset (1000 frames, transmit angle #3, 64 channels)
        ② Pre-computed Field II impulse-response library (64 elements × 36480 pixels)

── Stage 1: system definition (continue1) ──
c = 1540, f0 = 7.813 MHz, fs = 31.25 MHz, 64 elements, pitch = 0.2 mm
probe params ← getparam('L12-3v')                       # MUST toolbox
load measured impulse response (tv.mat, downsampled + rescaled)   # TX/RX characteristics
load measured transmit waveform (TW_L12_067DC_2cycles.mat)
define imaging grid: 128 (lateral) × 285 (depth) = 36480 pixels
             lateral ±32.73λ (≈ ±6.45 mm), depth 10.15λ–152.2λ (≈ 2–30 mm)

── Stage 2: observation matrix (continue2) ──
for frame f = 1..1000:
    take the RF segment of angle #3, channels 65:128 (len×64), NaN→0
    map to the 1200-point sub-band centered on f0 → unfold into a 76800×1 column
stack: Y = rf_all ∈ ℂ^(76800×1000)          # 76800 = 1200×64

── Stage 3: system matrix (continue3) ──
for element e = 1..64:
    load the impulse responses of this element for all 36480 pixels (36480 columns)
    same sub-band extraction → 1200×36480
vertically stack the 64 elements: R ∈ ℂ^(76800×36480)   # same observation space as Y

── Stage 4: TwIST sparse solve (TWIST_vera_L12_64eles) ──
τ  ← 0.02 · max|R'Y|
X₀ ← 0 (36480×1000)
X* ← TwIST(Y, R, τ;  Psi = soft2, Phi = weighted_L1, λ₁ = 1e-4,
           Monotone = 1, StopCriterion = 1, tol = 1e-2, Maxiter = 50)

── Stage 5: output (main script) ──
for f = 1..1000: img_f = reshape(|X*(:,f)|, 285, 128), linearly rescaled to [0,255]
write a 3D tif (285×128×1000) and save X*.mat
```

**Call chain**: main script → continue1 (parameters) / continue2 (RF→Y) / continue3 (impulse responses→R) → TWIST_vera (assembles the optimization problem) → TwIST.m (solver; uses soft2 as the proximal operator and weighted_L1 for objective evaluation) → main script post-processing/output.

---

## 3. Exact Objective Function

`TwIST.m` evaluates the objective `f = 0.5*resid'resid + tau*phi(x)` at every iteration; the problem actually solved by this code is:

$$\min_{X\in\mathbb{C}^{36480\times 1000}}\; F(X)=\tfrac{1}{2}\big\|Y-RX\big\|_F^2+\tau\,\Phi(X)$$

$$\Phi(X)=\underbrace{\big\|\sigma(X)\big\|_1}_{\text{temporal-std-map sparsity}}\;+\;\underbrace{\tfrac{1}{T}\big\|X(:)\big\|_1}_{\text{per-frame-averaged global L1}},\qquad \tau = 0.02\times\max\big|R'Y\big|$$

with the per-pixel temporal standard deviation (T = 1000 frames, unbiased N−1 normalization, complex deviations taken by magnitude):

$$\sigma_p=\sqrt{\tfrac{1}{T-1}\sum_{t=1}^{T}\big|X_{p,t}-\bar x_p\big|^2},\qquad \bar x_p=\tfrac{1}{T}\sum_{t=1}^{T}X_{p,t}$$

i.e., expanded:

$$\Phi(X)=\sum_{p=1}^{36480}\sigma_p\;+\;\frac{1}{T}\sum_{p=1}^{36480}\sum_{t=1}^{T}\big|X_{p,t}\big|$$

The first term enforces **sparsity of temporal fluctuation (vascular structure)**; the second enforces **spatio-temporal amplitude sparsity (microflow trajectory)**. Both are at per-pixel/per-element scale and are added directly.

(The λ₁ = 1e-4 passed via `'Lambda'` does **not** enter the objective; it is only used to compute the TwIST step-size parameters α and β — see Section 8.)

---

## 4. Parameter Dimensions

| Symbol / variable | Dimensions | Description |
|---|---|---|
| `Y` = `rf_all(:,1:1000)` = y | **76800 × 1000**, complex | Observation/RF matrix; 76800 = 1200 sub-band points × 64 receive elements; columns = slow-time frames |
| `R` = `simData_all` | **76800 × 36480**, complex | System matrix; rows live in the same observation space as Y, columns = pixels |
| `x_twist` = X* | **36480 × 1000**, complex | Optimization variable / reconstructed image sequence; 36480 = 285 depths × 128 lateral (pixels ordered as lateral blocks of 285 depths each) |
| `initarray` | 36480 × 1000 | All-zero initialization |
| `hRt(y)` = R'Y, `grad` | 36480 × 1000 | Gradient / back-projection |
| `resid` = Y−RX | 76800 × 1000 | Residual |
| `tau` | scalar 0.02 | Penalty factor; effective penalty = `tau·max(max(abs(R'Y)))` (per-pixel τ of size(x) is also supported but not used here) |
| `lambda1` | scalar 1e-4 | Affects only α, β; not in the objective |
| `bm_kk` | 285 × 128 | Single-frame magnitude image |
| `temp_bm_kk` | 285 × 128 × 1000, single | Output tif volume |
| `xstack` inside weighted_L1 | 285 × 128 × 1000 | `reshape(x,285,128,[])`, 3rd dim = frames |
| `x_mip2` = σ inside weighted_L1 | 285 × 128 | Temporal-fluctuation map after std over frames |
| `y2` = ‖σ‖₁, `y3` = ‖X(:)‖₁, `y3_scale` = y3/T | scalars | Regularizer intermediate quantities |

Physical grid: lateral ±6.45 mm, depth ≈ 2–30 mm; each frame's `RcvDataFrame` uses channels `65:128` (64 channels) and angle `3:3` (only one transmit angle).

---

## 5. Standard-Deviation Computation in the Regularizer

`02_optimization\weighted_L1.m`:

```matlab
xstack  = reshape(x, 285, 128, []);   % pixels × lateral × frames (slow time)
x_mip2  = std(xstack, 0, 3);          % per-pixel std along the 3rd dim (frame/slow-time axis)
y2      = norm(x_mip2, 1);            % ‖σ‖₁
```

- `std(xstack,0,3)`: **per-pixel temporal standard deviation along the 3rd dimension (frame/slow-time axis)**; weight flag 0 means normalization by **N−1 (=999)**, i.e., the unbiased estimator:

$$\sigma_p=\sqrt{\tfrac{1}{T-1}\sum_{t=1}^{T}\big|x_{p,t}-\bar x_p\big|^2}$$

- For complex data, MATLAB's `std` sums **magnitudes of deviations** from the complex mean x̄_p, so σ_p measures the amplitude fluctuation of each complex pixel time series.
- The resulting 285×128 σ map is essentially a **power-Doppler-like temporal-fluctuation map** (sensitive to flow/vasculature); `y2 = ‖σ‖₁` is the "temporal-fluctuation sparsity" constraint and constitutes the **first (active) term** of Φ.

---

## 6. Scaling of Each Regularization Term

`02_optimization\weighted_L1.m`:

| Term | Expression | Scaling | Status |
|---|---|---|---|
| **y2** | $\|\sigma(X)\|_1$ (L1 of the temporal-std map; vascular-structure sparsity) | none | **active** (first term) |
| y3 | $\|X(:)\|_1$ (spatio-temporal global L1; microflow-trajectory sparsity) | none | intermediate |
| **y3_scale** | $y_3/T$, T = `size(x,2)` = 1000 | **divided by the frame count** | **active** (second term) |
| **Φ = y** | $y_2 + y_3/T$ | direct sum of the two terms | final regularizer |

**Purpose of the scaling**: $\|X(:)\|_1$ grows linearly with the frame count T (more frames reconstructed → larger term), whereas $\|\sigma\|_1$ is essentially independent of T (the std of similar-amplitude time series does not grow with T). Dividing y3 by T brings both terms to a "per-frame scale", so the combined weighting no longer depends on the number of frames — i.e., a frame-count-independent balance between "structure sparsity" and "trajectory sparsity".

---

## 7. Use of the Penalty Parameter τ

**Value (`01_main_pipeline\TWIST_vera_L12_64eles.m`)**:

```matlab
tau = 0.02;                                        % penalty factor (named constant)
TwIST(y, hR, tau.*max(max(abs(hRt(y)))), ...)      % effective penalty parameter
```

$$\tau_{\text{eff}} = 0.02\times\max\big|R'Y\big|$$

i.e., **2%** of the matched-filter / back-projection peak: the factor 0.02 is extracted as a named constant and multiplied inline by the data-amplitude peak max|R'Y|, so the choice remains dimensionless and self-adaptive to data amplitude (scale the data by k and τ_eff scales by k). max|A'y| is exactly the LASSO zero-solution threshold — τ ≥ max|A'y| would yield the zero solution; 2% stays safely inside the non-trivial regime. (TwIST's internal zero-solution shortcut only activates for the default soft/L1 setup; here a custom Psi is passed, so that branch is skipped.)

**Three points of use within the algorithm**:

1. Objective function: `f = 0.5‖resid‖² + tau*Phi(x)` (both stopping criterion 1 and the monotonicity test are based on it);
2. Proximal-operator threshold: the IST step `x = Psi(xm1 + grad/max_svd, tau/max_svd)` — the threshold scales with the adaptive step as **τ/max_svd**;
3. Not in the gradient: `grad = AT(resid)` is independent of τ.

A per-pixel matrix τ (same size as x) is also supported; this code uses a scalar.

---

## 8. TwIST Implementation Details

### 8.1 Iteration structure (TwIST.m:481-529) — two-step iterative shrinkage/thresholding

```
init: xm1 = xm2 = X₀ = 0;  max_svd = 1;  Aty = A'y precomputed once
each outer iteration:
    grad ← AT(resid)                      # residual back-projection; not updated inside the inner loop
    inner (backtracking) loop:
        x ← Psi(xm1 + grad/max_svd, tau/max_svd)        # IST step
        if ≥2 IST steps done or already in TwIST phase:
            [sparse=1] mask = (x≠0); xm1,xm2 ← xm1.*mask, xm2.*mask   # zero the past where present is zero
            xm2 ← (α−β)·xm1 + (1−α)·xm2 + β·x          # two-step extrapolation
            f ← 0.5‖y−A xm2‖² + τΦ(xm2)
            if f increases and Monotone=1: fall back to IST (TwIST_iters = 0)
            else accept xm2, TwIST_iters++, break
        else (first two pure-IST warm-up steps):
            f ← 0.5‖y−A x‖² + τΦ(x)
            if f increases: max_svd ← 10·max_svd (shrink step, reset counters)   # handles ‖A‖²>1 divergence risk
            else break
    xm2 ← xm1;  xm1 ← x
    stopping: criterion = |f−f_prev|/f_prev (StopCriterion = 1)
          continue = (iter ≤ 50) and (criterion > 1e-2); at least 5 iterations
```

### 8.2 Step-size parameters (computed from λ₁ = 1e-4, λ_N = 1; TwIST.m:320-326)

- ρ₀ = (1−λ₁/λ_N)/(1+λ₁/λ_N) ≈ 0.99980
- **α = 2/(1+√(1−ρ₀²)) ≈ 1.9608**
- **β = 2α/(λ₁+λ_N) ≈ 3.9212**
- The two-step extrapolation coefficients are (α−β) ≈ −1.96, (1−α) ≈ −0.96, β ≈ 3.92 (typical underdamped extrapolation)

### 8.3 Further notes

- **Adaptive max_svd**: initialized to 1 (assuming ‖A‖₂² ≤ 1, which does not hold for the unnormalized R); ×10 back-off when monotonicity breaks during the IST phase; ×0.9 every 10000 TwIST iterations (slowly enlarging the step). The IST proximal threshold τ/max_svd moves in sync with the step 1/max_svd.
- **soft2 is the complex soft-threshold**: `y = max(|x|−T,0); y = y./(y+T).*x`, algebraically identical to `x·(1−T/|x|)₊`; it avoids division by zero and shrinks complex values by magnitude while preserving phase.
- **Division of labor between Phi and Psi**: `weighted_L1` (Phi) is used only to evaluate the objective (stopping + monotonicity test); the solution is actually shaped by `soft2` (Psi). Note that soft2 is the exact proximal operator of τ‖·‖₁ and does not strictly match the frame-coupled Φ — the shrinkage behavior of the iteration is equivalent to L1 soft-thresholding, while Φ (with its std term) influences stopping and monotone fallback through the objective value, thereby affecting the iteration progress of all 1000 frames (gradient and matrix products remain column-wise).
- **Debias = 0**: no CG debiasing phase; the output is the main-algorithm solution.
- **Verbose = 1**: prints objective/nonzero count at every iteration.
- The stopping criterion is the relative change of the **total objective over all 1000 frames** < 1%; convergence of individual low-amplitude frames is governed by the aggregate.

---

## 9. Usage

Add this folder (including all subfolders) to the MATLAB path and make sure the data files resolve:

```matlab
addpath(genpath('...\PER_L12_64eles_veraRF_recon_deps'));            % genpath includes all 5 subfolders automatically
addpath(genpath('K:\scui code copy ver-20250602\fieldii\pulse-echo impulse response\MPI'));  % tv.mat / TW_L12_067DC_2cycles.mat
run PER_L12_64eles_veraRF_recon.m
```

External data required at runtime (not functions; not distributed with this folder):

| Data | Purpose | Location |
|---|---|---|
| `tv.mat` | Measured impulse response (continue1) | `..\L12hydrophone\`, `..\MPI\` (resolved via the MATLAB path) |
| `TW_L12_067DC_2cycles.mat` | Measured transmit waveform (continue1) | `..\MPI\` |
| Verasonics RF dataset | Observation data (continue2) | `O:\OPTI\...`, `I:\OPTI\...` (absolute paths) |
| Simulated impulse-response matrices | System matrix (continue3) | `I:\OPTI\Normc_normalized_...` (absolute path) |

MATLAB built-in / toolbox functions (`rescale`, `downsample`, etc.) need no copying.

