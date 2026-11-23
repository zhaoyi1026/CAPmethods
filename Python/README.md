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

## What's new in 0.2.0

capcov 0.2.0 follows the R package CAPmethods 0.2.0 (HDCAP V9, CAP-CoC V4):

- **HDCAP / CAP: the objective.** The initial values of a direction are compared
  by the function the fit minimizes, `Σᵢ Tᵢ (xᵢ'β + sᵢ e^(−xᵢ'β))` at `γ'Hγ = 1`:
  `sᵢ = γ'Sᵢγ` without shrinkage (the negative log-likelihood), and with
  shrinkage the shrunk score of the β fit, `ρ₁γ'γ + ρ₂γ'Σᵢγ` (`Σᵢ` the
  per-subject `cov_ls`, `ρ₁, ρ₂` the fit's weights at β). 0.1.0 used the
  objective of the earlier R versions (shrinkage weights of another formula on
  the sample covariances, also for a fit without shrinkage). New:
  `capcov.obj_func(Y, X, gamma, beta, cov_shrinkage)`.
- **HDCAP / CAP: several initial values of β** — argument `beta0` of `cap_reg`,
  `cap_beta` and `cap_beta_boot`: a vector of length q or a `q × m` array (one
  initial value per column). The default is two, `0` and
  `capcov.cap_beta0_v5(q)`, R's `set.seed(100); c(0, rnorm(q − 1, 1, 1))` — the
  same numbers as in R, drawn without changing numpy's random numbers. Every
  initial value of γ runs from every initial value of β, and every β fit keeps
  the converged estimate with the lowest objective. 0.1.0 started from 0 only
  (`beta0=np.zeros(q)`).
- **HDCAP / CAP: initial values of γ** — `cap_reg(gamma0_mat=...)` (`p × k`); by
  default at most 10, as in R (0.1.0 used all `p + 6` random initial values).
- **CAP-CoC.** The objective that compares initial values uses the outcome
  covariances the fit iterates on (with outcome shrinkage `cov_sk_y` at the
  fitted variances, not the per-subject `cov_ls`); several initial values of β,
  `beta0` of `coc_reg`, `coc_coef`, `coc_coef_asmp` and `coc_coef_boot` (default
  `0` and `capcov.coc_beta0_v3(r)`, R's numbers; α starts at 0; 0.1.0: `0`);
  `coc_reg(gamma0_mat=..., theta0_mat=...)`.
- **All methods: failures.** A first direction that cannot be estimated raises a
  `RuntimeError` that names the step (0.1.0: a `TypeError` on `None`); a later
  direction that cannot be estimated ends the search, as R's `try()` does; an
  initial value whose fit or objective fails is left out.

Fits can differ from 0.1.0 where the new objective or initial values select
another candidate.

## Install

```bash
cd CAPmethods/Python          # in a clone of the repository
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

# initial values: beta0 (default 0 and capcov.cap_beta0_v5(q)), gamma0_mat (p x k)
fit0 = capcov.cap_reg(Y, X, nD=2, beta0=np.zeros(X.shape[1]))   # one initial value of beta
capcov.obj_func(Y, X, fit["gamma"][:, 0], fit["beta"][:, 0])     # objective of (gamma, beta)

# bootstrap inference for beta at a fixed direction
inf = capcov.cap_beta_boot(Y, X, fit["gamma"][:, 0], cov_shrinkage=False, sims=500)
inf["estimate"], inf["se"], inf["pvalue"], inf["lower"], inf["upper"]
```

`cap_reg` returns the directions; `cov_shrinkage` is auto-enabled when
`min_i T_i - 5 < p`. A first direction that cannot be estimated raises
`RuntimeError`; a later one ends the search.

## Usage (CAP-CoC)

```python
import numpy as np, capcov

# Y, X: lists of (T_i x q) / (T_i x p) arrays; W: (n x r) covariates incl. intercept
fit = capcov.coc_reg(Y, X, W, stop_crt="nD", nD=2)        # or stop_crt="DfD", DfD_thred=2
fit["gamma"], fit["theta"]   # q x k / p x k projection loadings (Y and X sides)
fit["alpha"], fit["beta"]    # k effects of log(theta'Sx theta) and r x k covariate effects
fit["score_y"], fit["score_x"]      # n x k projected variances
fit["DfD_y"], fit["DfD_x"]          # deviation-from-diagonality per direction

# inference for (alpha, beta) at a fixed (gamma, theta):
capcov.coc_coef_asmp(Y, X, W, gamma, theta)     # asymptotic
capcov.coc_coef_boot(Y, X, W, gamma, theta, sims=500)   # bootstrap

# initial values of beta (default 0 and capcov.coc_beta0_v3(r)); alpha starts at 0
fit = capcov.coc_reg(Y, X, W, stop_crt="nD", nD=2, beta0=np.zeros(W.shape[1]))
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

`tests/` has the R-free regression tests. The R cross-checks (2026-11-23, R 4.6.1)
ran against the installed R package CAPmethods 0.2.0 on identical data, with the
same initial values of γ (`gamma0_mat`, `ninitial` = its number of columns) and
the default initial values of β.

HDCAP / CAP vs R `capReg`:

- `cov_ls` linear shrinkage: **bit-identical** (max|Δ| ≈ 1e-15).
- R's random numbers: `_core.r_runif` / `r_rnorm` reproduce `set.seed(s); runif(n)`
  / `rnorm(n)` (runif bit-identical, rnorm to ≈ 2e-16), so `cap_beta0_v5` is R's
  initial value of β.
- `obj_func` (with and without shrinkage) and `cap_beta`: identical (relative
  Δ ≈ 1e-15; β within 4e-9 without shrinkage, where the two initial values reach
  the same β within the tolerance).
- **HDCAP** (`cov_shrinkage=True`, two directions; a simulation with p = 8 and
  the package example with p = 20): γ, β and scores **identical** (max|Δ| ≈ 1e-13).
- **CAP** (`cov_shrinkage=False`, two directions): direction 1 identical
  (max|Δ| ≈ 2e-10); direction 2 within 2e-3 in γ (cosine 0.999994) and 7e-5 in β.
  Before the second direction the data are deflated with an SVD whose basis for
  the removed direction is arbitrary, and R's LAPACK and numpy's choose different
  vectors; on R's deflated data the Python search returns R's direction exactly.

CAP-CoC vs R `COCReg`:

- `cov.ls` / `cov.sk.x` / `cov.sk.y` covariance kernels: **bit-identical** (max|Δ| ≈ 1e-15).
- **`COCReg.coef`** (α, β given γ, θ; with and without shrinkage), the objective
  with outcome shrinkage, and **`COCReg.coef.asmp`** (estimates, SE, p-values,
  CI): identical (max|Δ| ≈ 1e-15).
- **`COCReg`** on the package example with its recovering settings (`Hy = I`,
  `Hx = I`, `burn_in = 200`, nD = 1): γ identical (7e-16); θ, α, β within 6e-5
  (the convergence tolerance).
- On a simulation where CoC is weakly identified with the default `Hy`, `Hx`
  (the one initial value that finds the true direction has the highest
  objective), some initial values are numerically unstable: a relative change
  of 1e-12 in the start changes the direction they reach. R and Python reach
  different directions from these, so `COCReg` can return other directions; the
  stable initial values agree with R (≈ 1e-6).

The LCAP, CAP-mediation and CAP-clustering code of CAPmethods 0.2.0 gives the
results of the versions compared below (bit-identical in R), so these
comparisons carry over; capcov 0.2.0 adds their handling of failures.

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
                #         gamma_solve (generalized eigenproblem), diag_level (DfD),
                #         r_runif / r_rnorm (R's random numbers)
  hdcap.py      # HDCAP/CAP: cap_reg, cap_beta, cap_beta_boot, obj_func, cap_beta0_v5  [done]
  coc.py        # CAP-CoC: coc_reg, coc_coef, coc_coef_asmp, coc_coef_boot, coc_beta0_v3  [done]
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

The citations are the R package's. Access them programmatically via
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
- **CAP-HDcov** (`hdcov`): He, Y., Zou, C., & Zhao, Y. (2026). Covariance
  regression with high-dimensional predictors: An application to link brain
  structural and functional connectivity. *The Annals of Applied Statistics*,
  20(2), 1477–1495. <https://doi.org/10.1214/26-AOAS2157>
- **MCAP**, **CAP-clustering**: no publication yet.
