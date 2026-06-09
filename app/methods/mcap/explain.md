### What MCAP does

**Multilevel Covariate Assisted Principal (MCAP) regression** extends CAP to
data with a **hierarchical / multilevel** structure: units are nested in
**clusters** (e.g. age groups), and each cluster is allowed its **own** projection
direction. Where LCAP holds the projection $\gamma$ fixed across subjects, MCAP
lets it **vary across clusters** — each cluster $i$ has a direction $\gamma_i$ on
the unit sphere drawn from a **von Mises–Fisher** distribution centered on a
population direction $\gamma$ with concentration $\kappa$.

For cluster $i$, unit $j$, with covariance $\Sigma_{ij}$ of the $T_{ij}\times p$
data, the model (Equation (1) of the reference) is

$$\log\!\left(\gamma_i^\top \Sigma_{ij}\, \gamma_i\right)
   = \beta_0 + x_{1i(j)}^\top \beta_1 + x_{2ij}^\top(\beta_2 + \vartheta_i) + \varepsilon_i,
   \qquad \gamma_i \sim \mathrm{vMF}(\gamma, \kappa),$$

where $\beta_0$ is the intercept, $x_{1i(j)}$ are covariates with **fixed**
effects $\beta_1$, and $x_{2ij}$ are covariates with fixed effects $\beta_2$ and
**cluster-level random slopes** $\vartheta_i \sim N(0,\Omega)$;
$\varepsilon_i \sim N(0,\sigma^2)$ is a cluster **random intercept**. A larger
$\kappa$ means the cluster directions cluster tightly around $\gamma$; small
$\kappa$ means the covariance-associated subnetwork reorients markedly across
clusters.

### Inputs

- **Response** — a **nested list**: element $i$ is a list over cluster $i$'s units,
  and unit $j$ is a $T_{ij}\times p$ matrix. Upload as `.rds`/`.RData`, or a long
  CSV with columns `[cluster, unit, V1 … Vp]`.
- **Fixed-effect covariates** $x_1$ — a list whose $i$-th element is an
  $n_i \times q_1$ matrix (one row per unit), or a long CSV
  `[cluster, unit, cov…]`. These enter with **population (fixed)** effects $\beta_1$.
- **Random-effect covariates** $x_2$ — a list whose $i$-th element is an
  $n_i \times q_2$ matrix, or a long CSV `[cluster, unit, cov…]`. These carry the
  **cluster-level random slopes** $\vartheta_i$. No intercept column for either —
  the intercept $\beta_0$ is added internally.

### Key parameters

- **Number of directions** — choose either a **fixed number** $n_D$ or the **DfD
  criterion** (keep adding directions while the deviation-from-diagonality stays
  below a threshold). The leading direction is the most reliably recovered.
- **DfD threshold** — used when selecting by DfD.
- **# random initializations** — more restarts make the non-convex,
  directional-statistics optimization more robust (≥8 recommended for clean
  recovery; more is slower).

### Outputs

- $\gamma$ — the **population** projection loadings per direction.
- $\gamma_i$ — **cluster-specific** directions; the app reports each cluster's
  cosine alignment to the population $\gamma$ (the vMF spread).
- $\kappa$ — the von Mises–Fisher concentration per direction.
- $\beta$ — population fixed covariate effects on the log projected-variance.
- **Random-effect variances** — the random-intercept variance $\sigma^2$ and the
  random-slope variances (diagonal of $\Omega$) per direction.
- **Scores** — each cluster-unit's projected variance vs. its fitted value.
- **DfD** across dimensions.

### Built-in example

The example has 20 clusters, each with about 50 units and roughly 80 samples per
unit, $p = 5$, and **two** covariate-driven directions (basis directions 2 and 4).
Each cluster gets its own projection $\gamma_i$ drawn around the population
direction with concentration $\kappa = 10$; there are two fixed-effect covariates
and one random-slope covariate. With the default $n_D = 2$ the method recovers
**both** population directions (cosine $\approx 0.99$ and $0.94$) and the fixed
effects match the truth; this run takes about 2–3 minutes. Choosing $n_D = 1$
recovers the leading direction faster, and the DfD criterion selects the number of
directions automatically.
