### What LCAP does

**Longitudinal Covariate Assisted Principal (LCAP) regression** extends CAP to
data where each subject's covariance matrix is observed **repeatedly over time**
(multiple *visits*). It finds a single **time-invariant** projection $\gamma$ —
the same loading vector for every subject and visit — whose log projected
variance follows a **linear mixed model** in the (visit-level) covariates.

For subject $i$ at visit $v$ ($v = 1,\dots,V_i$), with covariance $\Sigma_{iv}$
of the $T_{iv} \times p$ data, LCAP assumes a projection $\gamma$ such that
**Equation (2.2)** of the reference holds:

$$\log\left(\gamma^\top \Sigma_{iv}\, \gamma\right)
   = \beta_0 + x_{1i}^\top \beta_1 + x_{2iv}^\top(\beta_2 + \nu_i) + u_i,$$

where $\beta_0$ is the intercept; $x_{1i}$ are **time-invariant** covariates with
fixed effects $\beta_1$; $x_{2iv}$ are **time-varying** covariates with fixed
effects $\beta_2$ and subject-level **random slopes** $\nu_i \sim N(0,\Omega)$;
and $u_i \sim N(0,\sigma^2)$ is a subject **random intercept**. Writing the random
intercept $\beta_{0i} = \beta_0 + u_i$ and random slopes $\beta_{2i} = \beta_2 + \nu_i$,
each subject gets its own baseline and time-varying-covariate sensitivities while
the projection $\gamma$ stays **constant over visits**. (A time-varying covariate
whose random-slope variance is zero behaves like a time-invariant one.)

Because $\hat{\Sigma}_{iv}$ can be noisy or singular when $T_{iv}$ is small
relative to $p$, LCAP can apply **linear (Ledoit–Wolf style) covariance
shrinkage** to each subject-visit covariance, exactly as in HDCAP.

> The sibling method **MCAP** (multilevel CAP) allows the projection $\gamma$ to
> **vary** across clusters; LCAP here keeps $\gamma$ invariant.

### Inputs

- **Response** — a **nested list**: element $i$ is itself a list over that
  subject's visits, and visit $v$ is a $T_{iv} \times p$ matrix. Upload as
  `.rds`/`.RData`, or as a long CSV with columns `[id, visit, V1 … Vp]`.
- **Covariates** — a list whose $i$-th element is an $n_{V_i} \times q$ matrix
  (one row per visit), or a long CSV with columns `[id, visit, cov1 … covq]`.
  An intercept ($\beta_0$) is added automatically. Covariates entered here are
  treated as **time-varying** $x_{2iv}$ — each gets a random slope $\nu_i$ — so
  they should vary within subject across visits; a covariate held constant within
  a subject behaves as a time-invariant $x_{1i}$ (its random-slope variance ≈ 0).

### Key parameters

- **Covariance shrinkage** — apply linear shrinkage to each $\hat{\Sigma}_{iv}$
  (automatically enabled when $\min_{i,v} T_{iv} - 5 < p$).
- **Number of directions chosen by** — a fixed number (`nD`) or the **DfD
  criterion**, which keeps adding directions while the deviation-from-diagonality
  stays below a threshold.
- **# random initializations** — more restarts make the non-convex optimization
  more robust.
- **Bootstrap replicates** — subject-level bootstrap resamples used for inference
  on the fixed effects $\beta$ (set to 0 to skip and show estimates only).

### Outputs

- $\gamma$ — the time-invariant projection loadings per direction.
- $\beta$ — fixed covariate effects on the log projected-variance, with
  **bootstrap** standard errors, test statistics, and confidence intervals.
- **Random-effect variances** — the random-intercept variance $\sigma^2$ and the
  random-slope variances (diagonal of $\Omega$) per direction.
- **Shrinkage weights** $(\rho_1,\rho_2)$ per direction (when shrinkage is used).
- **Scores** — each subject-visit's projected variance vs. its fitted value.
- **DfD** across dimensions.

### Built-in example

The example has $n = 100$ subjects, each seen at about 5 visits with $\approx 100$
samples per visit, $p = 20$ responses, a shared time-invariant projection, two
within-subject covariates, and a subject-level random intercept. **Two**
directions satisfy the CAP model — basis directions 2 and 4 — whose log-variance
depends on the covariates; the rest carry only the random baseline. With
covariance shrinkage on, the default $n_D = 2$ recovers **both** covariate-driven
directions (loading cosine $\approx 0.98$ and $0.94$) and the fixed-effect
coefficients match the truth.
