### What CAP-HDcov does

**Covariate Assisted Principal regression with High-Dimensional Covariates
(HCAP)** is the CAP variant for the case where the number of **covariates** $q$
is large (possibly $q > n$). As in CAP, a shared projection $\gamma$ summarizes
each subject's response covariance, and the projected log-variance is regressed
on the covariates:

$$\log\!\left(\gamma^\top \Sigma_i\, \gamma\right) = x_i^\top \beta,
   \qquad \beta \in \mathbb{R}^q \text{ sparse,}$$

but here $\beta$ is **high-dimensional and sparse** — only a few covariates drive
each covariance direction. Direction estimation uses a Gamma-family lasso
(`cv.glmnet`) at each iteration; **inference** on which covariates matter uses
**sample-splitting with multi-split aggregation** (SPARE / SSHDI) over $B$
resamples, giving a per-covariate **selection frequency** and an aggregated
sparse estimate with p-values.

### Inputs

- **Response** — a list of length $n$; element $i$ is a $T_i \times p$ matrix.
  Upload as `.rds`/`.RData`, or a long CSV (subject-id column + $p$ response
  columns).
- **Covariates** — an $n \times q$ matrix (one row per subject, $q$ possibly
  large), or a CSV whose first column is the subject id. An intercept is added
  automatically.

### Key parameters

- **Number of directions** — a fixed count (`nD`) or the DfD criterion. Each
  direction adds a full inference pass, so more directions = longer runtime.
- **# random initializations** — random restarts for the non-convex direction
  estimation (≥ 2).
- **Multi-split resamples** ($B$) — sample-split + aggregation resamples for the
  selection inference; more is stabler but slower.

### Outputs

- $\gamma$ — projection loadings per direction.
- **β selection inference** — per covariate: the aggregated sparse estimate, its
  SD, p-value, and **selection frequency** across the $B$ splits (covariates with
  frequency > 0.5 are flagged as selected).
- **Selected covariates** — the chosen covariate set per direction.
- On the built-in example, a **signal-recovery** table shows which true signal
  covariates were recovered, and a β truth-vs-estimate table.
- Plots: γ loadings, covariate **selection frequency**, and a **−log₁₀ p**
  significance plot with a Bonferroni line.

### Built-in example

The example is high-dimensional in the covariates: $n = 100$ subjects, $p = 5$
responses with $T_i = 100$ samples each, and $q = 200$ covariates. **Two**
covariate-driven directions are built in, and each is driven by only a small
**sparse** subset of the 200 covariates. The sparse estimator recovers the loading
directions, and the multi-split selection inference identifies the true signal
covariates of each direction. Runs in about 30–40 seconds.

> **Performance note.** Unlike the covariance-bound CAP methods, CAP-HDcov is
> **glmnet-bound** (the cost is the repeated lasso/GLM fits over the $q$
> covariates, not the covariance algebra). The built-in example runs in roughly
> half a minute; larger $q$, more directions, or more resamples take longer.
