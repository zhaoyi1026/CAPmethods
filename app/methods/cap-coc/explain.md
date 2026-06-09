### What CAP-CoC does

**Covariance-on-Covariance (CoC) regression** relates two sets of
covariance-matrix data. For each subject $i$ it summarizes a **predictor**
covariance $\Delta_i$ by its variance along a direction $\theta$, and an
**outcome** covariance $\Sigma_i$ by its variance along a direction $\gamma$, and
regresses one (log) projected variance on the other while adjusting for
covariates $W_i$:

$$\log\!\big(\gamma^\top \Sigma_i\, \gamma\big)
   = \alpha \,\log\!\big(\theta^\top \Delta_i\, \theta\big) + \beta^\top W_i.$$

Both projection directions $\gamma$ (outcome) and $\theta$ (predictor) are
**learned from the data**; $\alpha$ measures how strongly the predictor's
projected variance drives the outcome's, and $\beta$ are the covariate effects.
$\Sigma_i$ and $\Delta_i$ are estimated with linear (Ledoit–Wolf) shrinkage.
Multiple direction pairs can be extracted (each deflated from the previous).

### Inputs

- **Outcome** — a list of length $n$; element $i$ is a $T_{iy} \times q$ matrix
  (the repeated measurements whose $q \times q$ covariance is the outcome).
  Upload as `.rds`/`.RData`, or a long CSV (subject-id + $q$ outcome columns).
- **Predictor** — a list of length $n$; element $i$ is a $T_{ix} \times p$ matrix
  (covariance is the predictor). `.rds`/`.RData` or long CSV.
- **Covariates** — an $n \times r$ matrix $W$ (one row per subject), or a CSV
  whose first column is the subject id. An intercept is added automatically.

### Key parameters

- **Number of directions** — a fixed count (`nD`) or the DfD criterion. The
  **leading** direction pair is the most reliably recovered; later pairs are
  harder (weaker signal after deflation).
- **# random initializations** — restarts for the non-convex joint $\gamma/\theta$
  search.
- **Bootstrap replicates** — subject-level bootstrap for inference on $\alpha$ and
  $\beta$ (0 to skip).

### Outputs

- $\gamma$ / $\theta$ — outcome / predictor projection loadings per direction.
- $\alpha$, $\beta$ — the covariance-on-covariance coefficient and covariate
  effects, with **bootstrap** standard errors, test statistics, and confidence
  intervals.
- A **CoC-fit** scatter of $\log(\gamma^\top\Sigma\gamma)$ vs the fitted
  $\alpha\log(\theta^\top\Delta\theta)+\beta^\top W$.
- **DfD** (deviation-from-diagonality) for the outcome and predictor sides.

### Built-in example

The example has $p = 10$ predictor dimensions, $q = 5$ outcome dimensions, and 150
subjects with 150 samples each. There are **two** covariance-on-covariance pairs:
predictor directions drive outcome directions with regression coefficients
$\alpha = 3$ and $\alpha = 2$. The default $n_D = 1$ recovers the leading pair —
$\hat\alpha \approx 3$ with the $\gamma$ and $\theta$ directions recovered to
cosine $\approx 1$.
