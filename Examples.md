# Running the CAP methods — worked examples

Every method ships a built-in synthetic-data generator (`*_example()`, mirroring
the "built-in example" in the Shiny app) that returns data in exactly the shape
the wrapper expects, plus a `truth` element so you can check the estimates. The
pattern is always:

```r
library(CAPmethods)
d   <- <method>_example()      # generate example data (+ ground truth in d$truth)
fit <- <method>(...)           # run the method on it
```

| Method | Fit wrapper | Example | Demonstrates |
|--------|-------------|---------|--------------|
| HDCAP | `hdcap()` | `hdcap_example()` | covariate → covariance magnitude (one direction) |
| LCAP | `lcap()` | `lcap_example()` | longitudinal, fixed + random effects |
| MCAP | `mcap()` | `mcap_example()` | multilevel, cluster-varying loadings |
| CAP-CoC | `coc_d1()` / `coc()` | `coc_example()` | covariance-on-covariance regression |
| CAP-mediation | `capmediation()` | `capmediation_example()` | covariance mediator |
| HCAP | `hcap()` | `hcap_example()` | high-dimensional covariates (sparse loadings) |
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

**Data:** `X` is an `n × q` covariate matrix; `Y_list` is a length-`n` list of
`Tᵢ × p` response matrices.

```r
d   <- hdcap_example()                  # 80 subjects, p = 6, one covariate-driven direction
fit <- hdcap(d$Y_list, d$X, stop.crt = "nD", nD = 1, cov.shrinkage = FALSE)

fit$beta      # coefficient on the variance scale  → ~ c(0, 0.9)
fit$gamma     # estimated loading direction (p-vector)
```

Recovers the truth: **β̂ = (−0.01, 0.90)** vs the true `(0, 0.9)`, and the loading
**γ cosine ≈ 0.997**. The right panel shows each subject's variance score rising
with the covariate, exactly as `β > 0` implies.

![HDCAP](man/figures/hdcap.png)

---

## LCAP — longitudinal CAP (time-invariant projection)

A single time-invariant loading γ with a mixed model on the log-variance:
fixed effects β for time-varying covariates plus random intercept (σ²) and
random slopes (Ω). Visit `v` of subject `i` contributes a `Tᵢᵥ × p` matrix.

**Data:** `Y` is a nested list (subject → visit → `Tᵢᵥ × p`); `X` is a list of
`nVᵢ × q` covariate matrices (one row per visit).

```r
d   <- lcap_example()                   # 60 subjects × 6 visits, p = 6
fit <- lcap(d$Y, d$X, stop.crt = "nD", nD = 1, cov.shrinkage = FALSE, verbose = FALSE)

fit$beta          # fixed effects (Intercept, time, dose) → ~ c(0, -0.5, 0.6)
fit$beta0.sigma2  # random-intercept variance
fit$beta1.Omega   # random-slope covariance
```

Fixed effects recover the truth: **β̂ = (−0.01, −0.55, 0.63)** vs `(0, −0.5, 0.6)`;
**γ cosine ≈ 0.999**.

![LCAP](man/figures/lcap.png)

---

## MCAP — multilevel CAP (cluster-varying loadings)

Units nested in clusters; each cluster `i` has its **own** loading `γᵢ` drawn
around a population direction γ via a von Mises–Fisher distribution (concentration
κ). The log-variance has cluster random intercept/slopes.

**Data:** `Y` is a nested list (cluster → unit → `Tᵢⱼ × p`); `X2` is a list of
`nᵢ × q₂` random-effect covariate matrices (no intercept column).

```r
d   <- mcap_example()                   # 12 clusters × 15 units, p = 5
fit <- mcap(d$Y, X1 = NULL, X2 = d$X2, data.type = "Y",
            stop.crt = "DfD", DfD.thred = 2, method = "CAP",
            H.type = "CAvgCov", Omega.diag = TRUE, verbose = FALSE)

fit$gamma       # population loadings (p × nD)
fit$gamma.rnd   # per-cluster loadings (p × nD × m)
fit$kappa       # estimated vMF concentration
```

Population loading **γ cosine ≈ 0.98**. The right panel shows each cluster's
loading cosine to the population direction (their spread reflects κ̂).

![MCAP](man/figures/mcap.png)

---

## CAP-CoC — covariance-on-covariance regression

Regresses the **variance of one set of matrices on the variance of another**:
`log(γ′ Syᵢ γ) = α · log(θ′ Sxᵢ θ) + β′ Wᵢ`. It jointly finds the response
direction γ, the predictor direction θ, the association α, and covariate effects β.
`coc_d1()` estimates a single component; `coc()` runs the full multi-direction
selection.

**Data:** `Y` and `X` are length-`n` lists of `T × q` and `T × p` matrices; `W` is
an `n × r` covariate matrix.

```r
d   <- coc_example()                    # 80 subjects, X is p = 10, Y is q = 8
fit <- coc_d1(d$Y, d$X, d$W)

fit$alpha    # covariance-on-covariance association  → ~ 0.8
fit$gamma    # response (Y) loadings
fit$theta    # predictor (X) loadings
fit$beta     # covariate effects on Y-variance
```

Recovers **α̂ = 0.88** (true 0.8), response **γ cosine ≈ 0.99**, predictor
**θ cosine ≈ 0.88**.

![CAP-CoC](man/figures/coc.png)

---

## CAP-mediation — covariance mediator

Tests whether the **variance of a multivariate mediator** along a direction θ
carries an exposure effect to a scalar outcome: exposure → `log(θ′ Σᵢ θ)` → `Y`.

**Data:** `X` is an `n × q` design (exposure + covariates); `M` is a length-`n`
list of `Tᵢ × p` mediator matrices; `Y` is a length-`n` outcome vector.

```r
d   <- capmediation_example()           # 60 subjects, p = 10
fit <- capmediation(d$X, d$M, d$Y, stop.crt = "nD", nD = 1, verbose = FALSE)

fit$theta   # mediator loading direction
fit$coef    # path coefficients (exposure → mediator variance → outcome)
```

Mediator direction **θ cosine ≈ 0.99**. The right panel shows the fitted
mediator-variance → outcome path.

![CAP-mediation](man/figures/mediation.png)

---

## HCAP — high-dimensional-covariate CAP

Handles **many subject covariates** (`q` large) with sparse, `glmnet`-based
selection, returning `nD` response loadings and de-biased inference for the
covariate effects. This is the most expensive method (cross-validated lasso per
direction-iteration plus bootstrap inference) — use modest sizes when exploring.

**Data:** `X` is an `n × q` covariate matrix (here `q = 60`); `Y_list` is a
length-`n` list of `Tᵢ × p` matrices (`p = 5`).

```r
d   <- hcap_example(n = 80, p = 5, q = 60, Ti = 60)   # ~1 min
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
d   <- cappcl_example()                 # 80 subjects, p = 6, 2 covariance regimes
fit <- cappcl(d$Y, d$X, d$W, ncluster = 2, stop.crt = "nD", nD = 1)

fit$gamma   # shared loading direction
fit$class   # estimated cluster membership
fit$beta    # per-cluster variance-covariate coefficients
```

The shared loading is recovered (**γ cosine ≈ 0.99**); the right panel shows the
two covariance regimes capPCL separates (variance score by estimated cluster).

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
