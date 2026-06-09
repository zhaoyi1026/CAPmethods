### What HDCAP does

**High-Dimensional Covariate Assisted Principal (HDCAP) regression** extends CAP
to settings where the data dimension $p$ is large relative to the number of
samples per subject $T_i$ — the regime where each subject's sample covariance
$\hat{\Sigma}_{i}$ is noisy or even singular.

Like CAP, it finds projections $\gamma$ whose variance follows a log-linear model
in the covariates,

$$\log\left(\gamma^\top \Sigma_i\, \gamma\right) = x_i^\top \beta,$$

but it stabilizes estimation by applying **linear (Ledoit–Wolf style) shrinkage**
to each subject's covariance matrix,

$$\Sigma_i^{\text{shrunk}} = \rho_1 I + \rho_2 \hat{\Sigma}_{i},$$

with weights $(\rho_1,\rho_2)$ estimated from the data. This makes the method
well-behaved even when $T_i < p$.

### Inputs

- **Response** — a **list of length $n$** (one element per subject); the $i$-th
  element is a $T_i \times p$ matrix ($T_i$ = samples in subject $i$, $p$ = data
  dimension). Upload as `.rds`/`.RData`, or a long-format CSV.
- **Covariates** — an $n \times q$ matrix, one row per subject (or a CSV whose
  first column is the subject id). An intercept is added automatically.

### Key parameters

- **Covariance shrinkage** — apply linear shrinkage to each $\hat{\Sigma}_{i}$
  (recommended, and automatically enabled when $\min_i T_i - 5 < p$).
- **Number of directions chosen by** — a fixed number (`nD`) or the **DfD
  criterion**, which keeps adding directions while the deviation-from-diagonality
  stays below a threshold.
- **# random initializations** — more restarts make the non-convex optimization
  more robust.
- **Bootstrap replicates** — number of subject-level bootstrap resamples used for
  inference on $\beta$ (set to 0 to skip and show estimates only).

### Outputs

- $\gamma$ — projection loadings per direction.
- $\beta$ — covariate effects on the log projected-variance, with
  **bootstrap** standard errors, test statistics, and confidence intervals.
- **Shrinkage weights** $(\rho_1,\rho_2)$ per direction.
- **Scores** — each subject's projected variance (raw and shrinkage-adjusted).
- **DfD** across dimensions.

### Built-in example

The example has $n = 100$ subjects, $p = 20$ responses with $T_i = 100$ samples
each, a shared eigenbasis, and one binary covariate (`group`). **Two** directions
satisfy the CAP model — basis directions 2 and 4 — whose log-variance depends on
the covariate, with group effects $-1$ and $+1$; every other direction has a
covariate-free (random) log-variance. With covariance shrinkage on, the default
$n_D = 2$ recovers **both** covariate-driven directions (loading cosine
$\approx 0.99$, estimated group effects $\approx +1$ and $-1$). Increasing $n_D$
(or using the DfD criterion) also extracts the high-variance covariate-free
directions, whose estimated group effect $\approx 0$ correctly shows no covariate
association.
