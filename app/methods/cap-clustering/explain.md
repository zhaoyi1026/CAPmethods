### What CAP-clustering does

**Parsimonious Clustering of Covariance Matrices (CAPclust)** is model-based
clustering for subjects whose data are **covariance matrices**, sharing a common
CAP projection $\gamma$. Each subject's covariance is summarized by its projected
log-variance $\log(\gamma^\top \Sigma_i\, \gamma)$; within a cluster this follows a
log-linear model, and cluster membership depends on covariates:

$$\text{within cluster } k:\quad \log(\gamma^\top \Sigma_i\, \gamma) = x_i^\top \beta_k,$$
$$\text{membership:}\quad P(\text{cluster } k \mid w_i) \;\propto\; \exp(w_i^\top \alpha_k).$$

The projection $\gamma$, the cluster-specific variance coefficients $\beta_k$, and
the membership coefficients $\alpha_k$ are estimated jointly by **EM** — the E-step
gives posterior cluster probabilities, and the M-step updates $\gamma$ by a
generalized eigenproblem, $\beta_k$ by Newton steps, and $\alpha$ by a penalized
multinomial regression. Multiple **components** (directions) can be extracted,
each with its own clustering profile.

### Inputs

- **Covariance data** — a list of length $n$; element $i$ is a $T_i \times p$
  matrix whose covariance is clustered. Upload as `.rds`/`.RData`, or a long CSV
  (subject-id column + $p$ columns).
- **Variance-model covariates** $x_i$ — an $n \times q_1$ matrix for the
  within-cluster log-variance model.
- **Membership covariates** $w_i$ — an $n \times q_2$ matrix that drives cluster
  membership. Intercepts are added to both automatically.

### Key parameters

- **Number of clusters (K)** — the number of mixture components.
- **Number of components (nD)** — fixed (`nD`) or the DfD criterion. The leading
  component is the most reliably recovered; later components carry weaker signal.
- **# random initializations** — restarts for the non-convex γ/clustering EM.
- **Bootstrap replicates** — subject-level bootstrap for $\beta$ / $\alpha$
  inference (0 to skip).

### Outputs

- **Cluster sizes** and the recovered **cluster assignment** per component.
- $\beta_k$ — within-cluster variance coefficients (with bootstrap SE/p/CI), and a
  plot of $\beta$ by cluster with confidence intervals.
- $\alpha_k$ — membership (multinomial) coefficients (cluster 1 = reference).
- $\gamma$ — projection loadings; a box plot of the projected log-variance by
  estimated cluster (showing the separation).
- **DfD** across components.

### Built-in example

The example has 100 subjects, $p = 50$, $T_i = 100$ samples each, and $K = 2$
clusters, with two covariate-driven components — **D2** (well separated:
$\beta_1 = (1,1,-1)$, $\beta_2 = (-1,-1,1)$) and **D4** (clusters differ only in
covariate signs, harder). With the default $n_D = 1$ the method recovers the
leading component **D2**: the estimated $\gamma$ matches the truth (cosine
$\approx 1$) and the recovered clusters agree with the truth ($\approx 0.98$).
