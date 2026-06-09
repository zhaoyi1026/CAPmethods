# Running the CAP methods — worked examples

Every method ships a built-in synthetic-data generator (`*_example()`) that
returns data in exactly the shape the wrapper expects, plus a `truth` element so
you can check the estimates. The
pattern is always:

```r
library(CAPmethods)
d   <- <method>_example()      # generate example data (+ ground truth in d$truth)
fit <- <method>(...)           # run the method on it
```

| Method | Fit wrapper | Example | Demonstrates |
|--------|-------------|---------|--------------|
| HDCAP | `hdcap()` | `hdcap_example()` | covariate → covariance magnitude (two directions) |
| LCAP | `lcap()` | `lcap_example()` | longitudinal, fixed + random effects |
| MCAP | `mcap()` | `mcap_example()` | multilevel, cluster-varying loadings |
| CAP-CoC | `coc()` | `coc_example()` | covariance-on-covariance regression (manuscript sim) |
| CAP-mediation | `capmediation()` | `capmediation_example()` | covariance mediator (GMed sim) |
| HCAP | `hcap()` | `hcap_example()` | high-dimensional covariates (regularization + post-selection inference) |
| CAP-clustering | `cappcl()` | `cappcl_example()` | clustering covariance patterns |

All figures below are produced by the code shown (base graphics, no extra
packages). γ "loadings" are sign-invariant, so estimates are sign-aligned to the
truth before plotting; cosine similarity to the truth is reported in each title.

---

## HDCAP — high-dimensional CAP

Relates a subject covariate to the **magnitude of variance** along a latent
direction γ: `log(γ′ Σᵢ γ) = xᵢ′ β`. With `cov.shrinkage = FALSE` this is the
classical CAP model; set `TRUE` (auto when `Tᵢ − 5 < p`) for the shrinkage
estimator.

**Data:** the manuscript's HD-shrinkage simulation (`210309/eg`): `p = 20`,
`n = 100` subjects, `Tᵢ = 100`, one binary covariate (`group`), with **two**
covariate-driven directions (basis cols 2 and 4, group effects `−1` and `+1`).
`X` is `n × q`; `Y_list` is a length-`n` list of `Tᵢ × p` response matrices.

```r
d   <- hdcap_example()                  # p = 20, two covariate-driven directions
fit <- hdcap(d$Y_list, d$X, stop.crt = "nD", nD = 2, cov.shrinkage = TRUE)

fit$gamma     # two estimated loading directions (p × 2)
fit$beta      # group effect ~ +1 / -1 on the two directions
```

With shrinkage on, `nD = 2` recovers **both** covariate-driven directions
(**γ cosine ≈ 0.99 / 0.99**, **β̂(group) ≈ +1 / −1**). Increasing `nD` (or the DfD
criterion) additionally extracts the high-variance covariate-free directions,
whose `β̂(group) ≈ 0` correctly indicates no covariate effect.

![HDCAP](man/figures/hdcap.png)

---

## LCAP — longitudinal CAP (time-invariant projection)

A time-invariant loading γ with a mixed model on the log-variance: fixed effects
β plus random intercept (σ²) and random slopes (Ω). Visit `v` of subject `i`
contributes a `Tᵢᵥ × p` matrix. The example is the manuscript's `p20_q3`
simulation: `p = 20`, 100 subjects × ~5 visits, two within-subject covariates,
with **two** covariate-driven directions (basis cols 2 and 4).

**Data:** `Y` is a nested list (subject → visit → `Tᵢᵥ × p`); `X` is a list of
`nVᵢ × q` covariate matrices (one row per visit).

```r
d   <- lcap_example()                   # p = 20, two covariate-driven directions
fit <- lcap(d$Y, d$X, stop.crt = "nD", nD = 2, cov.shrinkage = TRUE, verbose = FALSE)

fit$gamma         # two estimated loading directions (p × 2)
fit$beta          # fixed effects (Intercept, x1, x2) per direction
fit$beta0.sigma2  # random-intercept variance
fit$beta1.Omega   # random-slope covariance
```

With shrinkage on, `nD = 2` recovers **both** covariate-driven directions
(**γ cosine ≈ 0.98 / 0.94**) and the fixed effects match the truth. `nD = 1`
recovers the leading direction faster; the DfD criterion selects the count.

![LCAP](man/figures/lcap.png)

---

## MCAP — multilevel CAP (cluster-varying loadings)

Units nested in clusters; each cluster `i` has its **own** loading `γᵢ` drawn
around a population direction γ via a von Mises–Fisher distribution (concentration
κ). The log-variance has cluster random intercept/slopes.

**Data:** the manuscript's γ-varying simulation (`p5_q4_2-1`, case 1): `p = 5`,
two covariate-driven directions, `m = 20` clusters. `Y` is a nested list
(cluster → unit → `Tᵢⱼ × p`); `X1` is a list of `nᵢ × q₁` **fixed**-effect
covariate matrices, `X2` a list of `nᵢ × q₂` **random**-slope covariate matrices
(no intercept columns). The number of directions is chosen by `nD` or `DfD`.

```r
d   <- mcap_example()                   # 20 clusters, 2 directions, X1 + X2
fit <- mcap(d$Y, X1 = d$X1, X2 = d$X2, data.type = "Y",
            stop.crt = "nD", nD = 2, ninitial = 8, method = "CAP",
            H.type = "CAvgCov", Omega.diag = TRUE, verbose = FALSE)

fit$gamma       # population loadings (p × nD)
fit$gamma.rnd   # per-cluster loadings (p × nD × m)
fit$kappa       # estimated vMF concentration per direction
```

Both directions recover: population loading **γ cosine ≈ 0.99 / 0.96**, and the
fixed-effect `β` match the truth. `nD = 1` recovers the leading direction faster;
the `DfD` criterion selects the number of directions automatically.

![MCAP](man/figures/mcap.png)

---

## CAP-CoC — covariance-on-covariance regression

Regresses the **variance of an outcome covariance on the variance of a predictor
covariance**: `log(γ′ Σᵧᵢ γ) = α · log(θ′ Σₓᵢ θ) + β′ Wᵢ`. It jointly finds the
outcome direction γ, the predictor direction θ, the association α, and covariate
effects β. `coc()` runs the multi-direction selection (used here with `nD = 1`);
`coc_d1()` estimates a single component directly.

**Data:** `X` (predictor) and `Y` (outcome) are length-`n` lists of `Tₓ × p` and
`Tᵧ × q` matrices; `W` is an `n × r` covariate matrix.

The built-in example is the **case-1 simulation of Zhao et al. (Biometrics 2025)**:
`p = 10` predictor / `q = 5` outcome, with two covariance-on-covariance pairs —
predictor directions {1, 3} drive outcome directions {2, 4} with `α = (3, 2)` and
`β = ((1,−1),(−1,1))`. The recovering fit uses **identity weights** `Hy = diag(q)`,
`Hx = diag(p)` (the default average-covariance weight collapses to a background
direction), `burn.in = 200`, and covariance shrinkage; `nD = 1` recovers the
leading pair (the deflated 2nd direction is unstable).

```r
d <- coc_example()                       # 150 subjects; X is p = 10, Y is q = 5
p <- ncol(d$X[[1]]); q <- ncol(d$Y[[1]])
fit <- coc(d$Y, d$X, d$W, stop.crt = "nD", nD = 1,
           Hy = diag(q), Hx = diag(p),
           cov.shrinkage.y = TRUE, cov.shrinkage.x = TRUE,
           burn.in = 200, ninitial = 5, seed = 100, verbose = FALSE)

fit$alpha    # covariance-on-covariance association  → ~ 3 (leading pair)
fit$gamma    # outcome (Y) loadings
fit$theta    # predictor (X) loadings
fit$beta     # covariate effects  → ~ c(1, -1)
```

Recovers the leading pair: **α̂ = 2.89** (true 3), outcome **γ cosine ≈ 1.0**,
predictor **θ cosine ≈ 1.0**, and **β̂ = (1.20, −0.96)** vs the true `(1, −1)`.

![CAP-CoC](man/figures/coc.png)

---

## CAP-mediation — covariance mediator

Tests whether a treatment acts on an outcome **through the variance of a
multivariate mediator** along a latent direction θ: treatment → `log(θ′ Σᵢ θ)`
(the **a-path**) → `Y` (the **b-path**), alongside a direct treatment effect.
`coef` reports `alpha` (a-path), `beta` (b-path), `gamma`, `IE` (indirect /
mediation effect) and `DE` (direct effect) per direction; `*_boot()` adds
bootstrap inference.

**Data:** `X` is an `n × q` design (treatment + optional covariates); `M` is a
length-`n` list of `Tᵢ × p` mediator matrices; `Y` is a length-`n` outcome vector.

The built-in example has 100 subjects, a mediator of dimension `p = 10` measured
with `Tᵢ = 150` samples each, and a binary treatment that shifts the mediator's
variance along one latent direction θ; the outcome depends on that log-variance.
The true indirect effect is `α·β = 0.8 × 0.7 = 0.56`. `H` defaults to the average
mediator covariance and `Y.remove = FALSE`.

```r
d   <- capmediation_example()           # 100 subjects, p = 10, binary treatment
fit <- capmediation(d$X, d$M, d$Y, stop.crt = "DfD", DfD.thred = 2,
                    Y.remove = FALSE, verbose = FALSE)

fit$theta   # mediator loading direction(s)
fit$coef    # alpha (a-path), beta (b-path), gamma, IE, DE  per direction
```

The leading mediating direction is recovered: **θ cosine ≈ 0.99**, and the
treatment → mediator-variance **a-path α̂ ≈ 0.8** (true 0.8) — the right panel
shows the mediator's variance score rising under treatment. The bootstrap
companion `capmediation_boot()` is the intended route for indirect-effect
inference.

![CAP-mediation](man/figures/mediation.png)

---

## HCAP — high-dimensional-covariate CAP

Handles **many subject covariates** (`q` large) with sparse, `glmnet`-based
selection, returning `nD` response loadings and de-biased inference for the
covariate effects. This is the most expensive method (cross-validated lasso per
direction-iteration plus post-selection inference) — use modest sizes when exploring.

**Data:** `X` is an `n × q` covariate matrix (here `q = 200`); `Y_list` is a
length-`n` list of `Tᵢ × p` matrices (`p = 5`). Two response directions are
covariate-driven, each by a small sparse subset of the covariates.

```r
d   <- hcap_example()                                 # p = 5, q = 200 (~30-40 s)
fit <- hcap(d$X, d$Y_list, stop.crt = "nD", nD = 2, B = 5)

fit$gamma_est   # estimated loadings (p × nD)
fit$inference   # per-direction de-biased covariate inference
```

HCAP returns the `nD` directions as a basis of the signal subspace (any rotation
within it is equivalent), so compare the **2-D subspace**, not individual columns:
canonical correlations to the truth are **1.00 and 0.999** (essentially exact).
The figure shows the estimated basis after Procrustes-aligning it to the truth.

![HCAP](man/figures/hcap.png)

---

## CAP-clustering — clustering covariance patterns

Clusters subjects whose covariance follows **different covariate models**, sharing
a common loading γ, with cluster membership informed by covariates `W`
(multinomial via `brglm2`). Each cluster gets its own variance-covariate
coefficients.

**Data:** `Y` is a length-`n` list of `Tᵢ × p` matrices; `X` is an `n × q₁`
variance-covariate matrix; `W` is an `n × q₂` membership-covariate matrix.

```r
d   <- cappcl_example()                 # 100 subjects, p = 50, K = 2 clusters
fit <- cappcl(d$Y, d$X, d$W, ncluster = 2, stop.crt = "nD", nD = 1)

fit$gamma   # shared loading direction
fit$class   # estimated cluster membership
fit$beta    # per-cluster variance-covariate coefficients
```

The leading component's shared loading is recovered (**γ cosine ≈ 1** to the true
direction D2); the right panel shows the two covariance regimes capPCL separates
(variance score by estimated cluster).

![CAP-clustering](man/figures/clustering.png)

---

### Notes

- **Bootstrap inference.** `hdcap_boot()`, `lcap_boot()`, `mcap_boot()`,
  `coc_boot()`, `capmediation_boot()`, and `cappcl_boot()` add bootstrap standard
  errors / CIs for the coefficients; pass `sims =` to set the number of replicates.
- **Direction selection.** Most fits accept `stop.crt = "nD"` (fix the number of
  directions via `nD`) or `stop.crt = "DfD"` (data-driven via a `DfD.thred`).
- **Any internal function** of a method is reachable with
  `cap_internal("<method>", "<fn>")` — e.g. `cap_internal("coc", "COCReg.D1")`.
