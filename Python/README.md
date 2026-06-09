# capcov — CAP family of covariance-regression methods (Python)

A Python package porting the **Covariate Assisted Principal (CAP)** family of
covariance-matrix regression methods from their R reference implementations.

> **CAP is not a separate module.** As requested, plain CAP is obtained from
> HDCAP with shrinkage off: `capcov.cap_reg(Y, X, cov_shrinkage=False)`. HDCAP
> with `cov_shrinkage=True` adds linear (Ledoit–Wolf) covariance shrinkage.

## Status

| Module | Method | Status |
|--------|--------|--------|
| `hdcap` | **HDCAP / CAP** | ✅ implemented, **verified vs R** |
| `coc` | **CAP-CoC** (covariance-on-covariance) | ✅ implemented, **verified vs R** |
| `lcap` | **LCAP** (γ-invariant longitudinal) | ✅ implemented, **verified vs R** |
| `mcap` | MCAP (γ-varying multilevel, vMF) | ⬜ planned (stub) |
| `mediation` | **CAP-mediation** (graph mediator) | ✅ implemented (approx.), **verified vs R** |
| `clustering` | **CAP-clustering** (PCL) | ✅ implemented (approx.), truth-verified |
| `hdcov` | CAP-HDcov (HCAP, high-dim covariates) | ⬜ planned (stub) |

The package is installable now with HDCAP/CAP, CAP-CoC, LCAP, CAP-clustering, and
CAP-mediation working; the remaining modules (MCAP, CAP-HDcov) import but raise
`NotImplementedError` with a pointer to the R version until ported (see roadmap).
HDCAP/CAP-CoC/LCAP are bit/optimizer-verified against R; CAP-clustering and
CAP-mediation are approximate ports (they substitute a self-contained solver for
an R package: `brglm2::brmultinom` and `nlme::lme` respectively) verified to agree
with R / recover the truth on the identifiable quantities.

## Install

```bash
cd CAP-Python/260605
pip install -e .            # core (numpy, scipy) -> HDCAP/CAP
pip install -e ".[full]"    # + scikit-learn / statsmodels / pandas (for the planned methods)
```

Requires Python ≥ 3.7, numpy ≥ 1.16, scipy ≥ 1.2.

## Built-in examples

Each implemented method has a self-contained example-data generator in
`capcov.examples` that returns data in the shape the wrapper expects plus a
`truth` entry (data-generating parameters), so estimates can be checked. These
mirror the R `CAPmethods` `*_example()` generators (manuscript / Shiny-demo
settings, kept in sync): HDCAP (p=20, two covariate-driven directions), LCAP
(p=20, two within-subject covariates), CoC (case-1 of Zhao & Zhao, *Biometrics*
2025), mediation (single-treatment sim of Xu & Zhao, *Biostatistics* 2025, with
IE = α·β = 0.8·0.7 = 0.56), and clustering (p=50, two components D2/D4).

```python
import capcov
from capcov import examples

d = examples.hdcap_example()                       # X, Y, truth
fit = capcov.cap_reg(d["Y"], d["X"], stop_crt="nD", nD=1, cov_shrinkage=False)
abs(fit["gamma"][:, 0] @ d["truth"]["gamma"])      # ~1.0 (direction recovered)

# others: examples.lcap_example(), coc_example(), mediation_example(),
#         clustering_example()  — each pairs with its wrapper (see Usage below).
```

## Usage (HDCAP / CAP)

```python
import numpy as np, capcov

# Y: list of (T_i x p) arrays (one per subject); X: (n x q) covariates incl. intercept
fit = capcov.cap_reg(Y, X, stop_crt="nD", nD=2, cov_shrinkage=False)  # plain CAP
fit = capcov.cap_reg(Y, X, stop_crt="DfD", DfD_thred=5, cov_shrinkage=True)  # HDCAP

fit["gamma"]          # p x nD projection loadings
fit["beta"]           # q x nD covariate effects on log projected-variance
fit["score"]          # n x nD subject projected variances
fit["DfD"]            # {'avg_level', 'sub_level'}

# bootstrap inference for beta at a fixed direction
inf = capcov.cap_beta_boot(Y, X, fit["gamma"][:, 0], cov_shrinkage=False, sims=500)
inf["estimate"], inf["se"], inf["pvalue"], inf["lower"], inf["upper"]
```

`cap_reg` returns the directions; `cov_shrinkage` is auto-enabled when
`min_i T_i - 5 < p`.

## Usage (CAP-CoC)

```python
import capcov

# Y, X: lists of (T_i x q) / (T_i x p) arrays; W: (n x r) covariates incl. intercept
fit = capcov.coc_reg(Y, X, W, stop_crt="nD", nD=2)        # or stop_crt="DfD", DfD_thred=2
fit["gamma"], fit["theta"]   # q x k / p x k projection loadings (Y and X sides)
fit["alpha"], fit["beta"]    # k effects of log(theta'Sx theta) and r x k covariate effects
fit["score_y"], fit["score_x"]      # n x k projected variances
fit["DfD_y"], fit["DfD_x"]          # deviation-from-diagonality per direction

# inference for (alpha, beta) at a fixed (gamma, theta):
capcov.coc_coef_asmp(Y, X, W, gamma, theta)     # asymptotic
capcov.coc_coef_boot(Y, X, W, gamma, theta, sims=500)   # bootstrap
```

Models `log(gamma' Sy_i gamma) = alpha * log(theta' Sx_i theta) + W_i' beta`;
both `cov_shrinkage_y`/`cov_shrinkage_x` auto-enable when `min_i T_i - 5 < dim`.

## Usage (LCAP, longitudinal)

```python
from capcov import lcap

# Y[i][v] is a (T_iv x p) matrix (subject i, visit v); X[i] is (nV_i x q) design
fit = lcap.cap_reg(Y, X, stop_crt="nD", nD=1, cov_shrinkage=False)
fit["gamma"]          # p x k shared projection loadings
fit["beta"]           # q x k fixed effects (intercept + slopes) of log proj-variance
fit["beta0_random"], fit["beta1_random"]   # per-subject random intercept / slopes

est = lcap.cap_beta(Y, X, fit["gamma"][:, 0], cov_shrinkage=False)  # random-effects fit
inf = lcap.cap_beta_boot(Y, X, fit["gamma"][:, 0], cov_shrinkage=False, sims=200)
```

Subject-level random-effects model
`log(gamma' Sigma_iv gamma) = x_iv' b_i`, `b_i ~ N(beta, diag(sigma2, Omega))`
(diagonal random-effect covariance, matching the web app).

## Usage (CAP-mediation)

```python
from capcov import mediation   # or capcov.cap_med

# X: (n x q) exposure+covariates; M: list of (T_i x p) mediator data; Y: (n,) outcome
fit = mediation.cap_med(X, M, Y, stop_crt="nD", nD=1)
fit["theta"]   # p x k mediator projection
fit["alpha"]   # exposure -> mediator (path a)
fit["beta"]    # mediator -> outcome (path b)
fit["gamma"]   # direct effect
fit["IE"]      # indirect (mediation) effect = alpha_x * beta
```

## Usage (CAP-clustering)

```python
from capcov import clustering

# Y[i] is a (T_i x p) matrix; X is (n x q1) variance covariates (incl. intercept);
# W is (n x q2) membership covariates (incl. intercept)
fit = clustering.cap_pcl(Y, X, W, ncluster=2, stop_crt="nD", nD=1)
fit["gamma"]      # p x k shared projection loadings
fit["beta"]       # list of (q1 x K) per-cluster variance effects
fit["class"]      # n x k cluster assignment per direction

cf = clustering.cap_pcl_coef(Y, X, W, fit["gamma"][:, 0], ncluster=2)  # at a fixed gamma
cf["class"], cf["logLik"]
```

## Verification

Validated against the R `capReg` on identical simulated data
(`tests/` has the R-free regression tests; the R cross-check script lives in the
project history):

- `cov_ls` linear shrinkage: **bit-identical** (max|Δ| ≈ 1e-15).
- **CAP** (`cov_shrinkage=False`): γ cosine = 1.000, β within ~1e-10 of R.
- **HDCAP** (`cov_shrinkage=True`): γ cosine ≈ 0.99999, β within ~1e-3 (the
  shrinkage-weight iteration path differs slightly between the RNG-seeded R init
  and the Python init; both converge to the same estimate and recover the truth).

CAP-CoC validated against R `COCReg` on identical simulated data:

- `cov.ls` / `cov.sk.x` / `cov.sk.y` covariance kernels: **bit-identical** (max|Δ| ≈ 1e-15).
- **`COCReg.coef`** (alpha/beta given gamma, theta): match within ~1e-7 (no
  shrinkage) and ~1e-6 (shrinkage).
- **`COCReg`** (`nD=2`): γ cosine = `[1.000, 0.9998]`, θ cosine = `[1.000, 1.000]`,
  α/β direction-1 essentially exact (β max|Δ| ≈ 1.6e-5), DfD within ~1e-3 (the
  small direction-2 gap is the shrinkage covariance-update path, as for HDCAP).
- **`COCReg.coef.asmp`** asymptotic inference (α, β SE / statistic / p-value /
  CI): identical to 4 decimals.

LCAP validated against R `capReg` on identical simulated data:

- **`cap_beta`** (random-effects fit at a fixed gamma), both no-shrinkage and
  shrinkage: β, σ², Ω all **bit-identical** (max|Δ| ≈ 1e-16).
- **`cap_D1`** (one direction at a fixed γ₀), no shrinkage: γ cosine = 1.000,
  β bit-identical; shrinkage: γ cosine ≈ 0.9996 (the shrinkage direction path
  recomputes the constant-shrinkage covariance each iteration, a slightly
  different path that converges to the same direction).
- **`capReg`** (`nD=1`): recovered γ cosine = 0.99998 vs R's direction 1.
  (Multi-start γ₀ uses numpy's RNG, so individual starts differ from R's, but the
  start-robust objective minimum agrees.)

CAP-clustering validated against R `capPCL`:

- Covariance kernels (`smat` second-moment, projected `score`, weighted `accum`,
  `eigen.solve`) and `diag.level`: **bit-identical** (max|Δ| = 0).
- The EM has many local optima; R's `capPCL_coef` uses a single `set.seed(100)`
  start, so it can land in a worse optimum (on the benchmark data Python's
  multi-restart fit reaches a strictly higher log-likelihood, −9272 vs −10797).
- On well-separated data the multi-restart EM **recovers the truth** (cluster
  agreement = 1.000, β within 0.02); `capPCL_D1` recovers the true γ (cos = 0.985).
- The membership model approximates `brglm2::brmultinom` (Firth multinomial) with
  a self-contained weighted Firth logistic, so membership *coefficients* may
  differ slightly from R; assignments/β/γ agree. This is a statistically-
  equivalent port, not bit-identical (R's EM is single-start + brmultinom).

CAP-mediation validated against R `CAPMediation_D1`:

- Covariance kernels (`med_cov`, projected `score`, weighted `accum`,
  `eigen.solve`): bit-identical.
- R fits the M-model with `nlme::lme(score ~ X, random=~1|ID)` on **singleton
  groups**, so the random-intercept/residual variance split is non-identifiable
  and R's BLUP is optimizer-arbitrary. This port keeps the identifiable fixed
  effects (GLS = OLS) and a deterministic BLUP (`blup_shrink`, default 1).
- The identifiable outputs match R closely: θ cosine = **0.9998**, α matches to
  ~3 dp, β within ~0.004, IE within ~0.003.
- The multi-start direction selection uses the objective at the **H-unscaled**
  theta (as R does), which penalizes high-variance background directions so the
  mediating direction is recovered rather than the dominant-variance one. On the
  bundled `mediation_example` (single-treatment sim) this recovers θ (cos ≈ 0.996)
  and the a-path (α̂ ≈ 0.8); the b-path/IE are mildly attenuated at finite Tᵢ.

Run the test suite:

```bash
pip install -e ".[test]" && pytest -q
```

## Architecture

```
capcov/
  _core.py      # shared: cov_ls shrinkage, sigma_cube, scores, accum,
                #         gamma_solve (generalized eigenproblem), diag_level (DfD)
  hdcap.py      # HDCAP/CAP: cap_reg, cap_beta, cap_beta_boot   [done]
  coc.py        # CAP-CoC: coc_reg, coc_coef, coc_coef_asmp, coc_coef_boot  [done]
  lcap.py       # LCAP: cap_reg, cap_beta, cap_beta_boot, diag_level  [done]
  clustering.py # CAP-clustering: cap_pcl, cap_pcl_coef, cap_pcl_coef_boot  [done]
  mediation.py  # CAP-mediation: cap_med, cap_med_coef, cap_med_d1  [done]
  examples.py   # built-in example-data generators (one per implemented method)
  references.py # per-method citations + DOIs; capcov.references()
  mcap.py       # MCAP (gamma-varying)     [stub]
  hdcov.py      # CAP-HDcov (HCAP)  [stub]
```

All methods share `_core`: the per-subject covariance, projected scores
`v'Σv`, weighted accumulation `Σ wᵢ Σᵢ`, and the generalized smallest-eigenvector
solve `gamma_solve(A, H)` (scipy `eigh(A, H)`).

## Porting roadmap (remaining methods)

Each will be ported from its R `Vn` source and cross-checked numerically. Notes
on Python dependencies for the non-linear-algebra pieces:

- **CAP-CoC** — ✅ done (pure numpy/scipy: two projections, alternating WLS +
  generalized eigen solve, linear/shrinkage covariance kernels).
- **LCAP** — ✅ done (subject-level random-effects EM recursion + generalized
  eigen direction update; numpy/scipy only).
- **MCAP** — γ-varying random-effects recursion; needs a von Mises–Fisher
  sampler/density (implement directly or via scipy).
- **CAP-mediation** — ✅ done (self-contained random-intercept fit; statsmodels
  intentionally avoided — not needed and segfaults in some envs).
- **CAP-clustering** — ✅ done (EM; membership via a self-contained weighted Firth
  logistic approximating `brglm2::brmultinom`; multi-restart to escape local optima).
- **CAP-HDcov** — Gamma-family lasso (`cv.glmnet`) + multi-split inference; the
  hard one (no drop-in Gamma-lasso in Python). Options: `glmnet`/`celer`-style
  solver or a custom Gamma-IRLS + coordinate descent.

The R implementations remain the reference; this package is being brought up to
parity method-by-method.

## References

Source of truth: `setting.md` in the project root. Access programmatically via
`capcov.references()` (all methods) or `capcov.references("coc")` (one method);
the data live in `capcov.REFERENCES`.

- **HDCAP / CAP** (`hdcap`, `cap_reg`):
  - Zhao, Y., Caffo, B., Luo, X., & ADNI (2021). Principal regression for high
    dimensional covariance matrices. *Electronic Journal of Statistics*, 15(2),
    4192. <https://doi.org/10.1214/21-EJS1887>
  - Zhao, Y., Wang, B., Mostofsky, S. H., Caffo, B. S., & Luo, X. (2021).
    Covariate assisted principal regression for covariance matrix outcomes.
    *Biostatistics*, 22(3), 629–645. <https://doi.org/10.1093/biostatistics/kxz057>
    (the classical CAP that HDCAP subsumes)
- **CAP-CoC** (`coc`): Zhao, Y., & Zhao, Y. (2025). Covariance-on-covariance
  regression. *Biometrics*, 81(3), ujaf097. <https://doi.org/10.1093/biomtc/ujaf097>
- **LCAP** (`lcap`): Zhao, Y., Caffo, B. S., & Luo, X. (2024). Longitudinal
  regression of covariance matrix outcomes. *Biostatistics*, 25(2), 385–401.
  <https://doi.org/10.1093/biostatistics/kxac045>
- **CAP-mediation** (`mediation`): Xu, Y., & Zhao, Y. (2025). Mediation analysis
  with graph mediator. *Biostatistics*, 26(1), kxaf004.
  <https://doi.org/10.1093/biostatistics/kxaf004>
- **MCAP**, **CAP-HDcov**, **CAP-clustering**: no publication yet.
