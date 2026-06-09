### What CAP does

**Covariate Assisted Principal (CAP) regression** finds linear projections of a
multivariate signal whose **variance** is associated with subject-level
covariates. It is a *covariance regression* method: instead of modeling the mean
of the responses, it models how the covariance structure of the responses
changes with covariates.

For each subject *i* we observe a $T_i \times p$ matrix of responses
$Y_i$ (e.g. $T_i$ time points over $p$ brain regions) and a vector of
covariates $x_i$. CAP estimates a direction $\gamma \in \mathbb{R}^p$ and
coefficients $\beta$ such that the projected variance follows a log-linear model:

$$\log\left(\gamma^\top \Sigma_i\, \gamma\right) = x_i^\top \beta,$$

where $\Sigma_i$ is subject *i*'s covariance matrix. The **score**
$\gamma^\top \Sigma_i \gamma$ is the variance of the data projected onto
$\gamma$; $\beta$ tells you which covariates drive that variance up or down.

### Inputs

- **Response** — a **list of length $n$** (one element per subject). The $i$-th
  element is a $T_i \times p$ matrix, where $T_i$ is the number of samples in
  subject $i$ and $p$ is the data dimension (e.g. the number of brain regions in
  an fMRI study). Upload it as an `.rds`/`.RData` file. *(A long-format CSV —
  subject-id column plus $p$ response columns — is also accepted and converted
  to this list internally.)*
- **Covariates** — an $n \times q$ matrix with one row per subject, in the same
  order as the response list (or a CSV whose first column is the subject id). An
  intercept (baseline log-variance) is added automatically.

### Key parameters

- **Number of directions chosen by** — either a **fixed number (`nD`)** of
  projections, or the **DfD criterion**, which keeps adding directions while the
  *deviation-from-diagonality* (DfD) stays below a threshold (e.g. DfD ≤ 2). The
  DfD measures how far the projected subject covariances are from being jointly
  diagonal; the search stops once an added direction would exceed the threshold.
- **Variant** — `CAP` (the optimization-based estimator) or `CAP-C` (common
  principal components flavor).
- **Orthogonal constraint (`OC`)** — enforce strict orthogonality between
  successive directions.
- **# random initializations (`ninitial`)** — more initializations make the
  non-convex optimization more robust at the cost of speed.

### Outputs

- $\gamma$ — the estimated projection loadings (per direction).
- $\beta$ — covariate effects on the log projected-variance, with asymptotic
  standard errors, test statistics and confidence intervals.
- **Scores** — each subject's projected variance, which you can relate back to
  the covariates.

### Built-in example

The example has $n = 80$ subjects, $p = 6$ responses with $T_i = 120$ observations
each, and two covariates: a binary `group` and a continuous `age`. **Two**
covariate-driven directions are built in:

- direction 1 — `group` effect $0.9$, `age` effect $-0.5$;
- direction 2 — `group` effect $-0.7$, `age` effect $0.4$.

With the default $n_D = 2$ (or the DfD criterion) the method recovers both: the
estimated loadings $\gamma$ match the truth (cosine $\approx 1$) and the estimated
coefficients match the values above.
