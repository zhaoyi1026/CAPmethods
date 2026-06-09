# Using `capcov`

A practical guide to running the Covariate-Assisted Principal (CAP) family of
covariance-regression methods in Python. For installation see the **Python**
section of the repository [README](../README.md); for the full method status and
verification notes see [`README.md`](README.md) in this folder.

Implemented methods (all verified against the R reference): **HDCAP/CAP**,
**CAP-CoC**, **LCAP**, **CAP-mediation**, **CAP-clustering**. (`mcap` and `hdcov`
import but raise `NotImplementedError` until ported.)

## The pattern

Every implemented method ships a built-in example-data generator in
`capcov.examples` that returns data in exactly the shape the function expects,
plus a `truth` entry (the data-generating parameters) so you can check estimates:

```python
import capcov
from capcov import examples

d   = examples.hdcap_example()        # dict: X, Y, truth
fit = capcov.cap_reg(d["Y"], d["X"], stop_crt="nD", nD=2, cov_shrinkage=False)
```

On your own data, pass arrays in the shapes documented per method below.

## HDCAP / CAP — `capcov.cap_reg`

Plain CAP is HDCAP with shrinkage off (`cov_shrinkage=False`); `cov_shrinkage=True`
adds linear (Ledoit–Wolf) shrinkage and auto-enables when `min_i T_i − 5 < p`.

```python
import capcov

# Y: list of (T_i x p) arrays, one per subject;  X: (n x q) covariates (incl. intercept)
fit = capcov.cap_reg(Y, X, stop_crt="nD", nD=2, cov_shrinkage=False)        # plain CAP
fit = capcov.cap_reg(Y, X, stop_crt="DfD", DfD_thred=5, cov_shrinkage=True) # HDCAP

fit["gamma"]   # p x nD  projection loadings
fit["beta"]    # q x nD  covariate effects on the log projected-variance
fit["score"]   # n x nD  per-subject projected variances
fit["DfD"]     # {'avg_level', 'sub_level'}

# bootstrap inference for beta at one fitted direction
inf = capcov.cap_beta_boot(Y, X, fit["gamma"][:, 0], cov_shrinkage=False, sims=500)
inf["estimate"], inf["se"], inf["pvalue"], inf["lower"], inf["upper"]
```

## CAP-CoC — `capcov.coc_reg`

Regresses the variance of an outcome covariance on the variance of a predictor
covariance: `log(γ' Sy_i γ) = α · log(θ' Sx_i θ) + W_i' β`.

```python
import numpy as np, capcov

# Y, X: lists of (T_i x q) / (T_i x p) arrays;  W: (n x r) covariates (incl. intercept)
fit = capcov.coc_reg(Y, X, W, stop_crt="nD", nD=2)        # or stop_crt="DfD", DfD_thred=2

fit["gamma"], fit["theta"]      # q x k / p x k  outcome- and predictor-side loadings
fit["alpha"], fit["beta"]       # k assoc. effects, r x k covariate effects
fit["score_y"], fit["score_x"]  # n x k projected variances

# inference for (alpha, beta) at fixed (gamma, theta)
capcov.coc_coef_asmp(Y, X, W, gamma, theta)           # asymptotic
capcov.coc_coef_boot(Y, X, W, gamma, theta, sims=500) # bootstrap
```

> On data like the bundled `coc_example` (two CoC pairs), pass identity weights
> `Hy=np.eye(q)`, `Hx=np.eye(p)` and a burn-in so the search finds the
> covariate-driven directions rather than a high-variance background one (the
> default average-covariance weight can collapse to background) — see the
> `coc_reg` signature for the `Hy`/`Hx`/`burn_in` arguments.

## LCAP — `capcov.lcap`

Longitudinal CAP with a time-invariant projection and a subject-level
random-effects model on the log projected-variance.

```python
from capcov import lcap

# Y[i][v] is a (T_iv x p) matrix (subject i, visit v);  X[i] is (nV_i x q) design
fit = lcap.cap_reg(Y, X, stop_crt="nD", nD=1, cov_shrinkage=False)

fit["gamma"]                              # p x k shared projection loadings
fit["beta"]                               # q x k fixed effects
fit["beta0_random"], fit["beta1_random"]  # per-subject random intercept / slopes

# random-effects fit and bootstrap at a fixed direction
est = lcap.cap_beta(Y, X, fit["gamma"][:, 0], cov_shrinkage=False)
inf = lcap.cap_beta_boot(Y, X, fit["gamma"][:, 0], cov_shrinkage=False, sims=200)
```

## CAP-mediation — `capcov.mediation`

Tests whether a treatment acts on an outcome through the variance of a
multivariate mediator along a latent direction.

```python
from capcov import mediation   # mediation.cap_med is also capcov.cap_med

# X: (n x q) exposure (+ covariates);  M: list of (T_i x p) mediator arrays;  Y: (n,) outcome
fit = mediation.cap_med(X, M, Y, stop_crt="nD", nD=1)

fit["theta"]   # p x k mediator projection
fit["alpha"]   # exposure -> mediator variance  (path a)
fit["beta"]    # mediator variance -> outcome   (path b)
fit["gamma"]   # direct effect
fit["IE"]      # indirect (mediation) effect
```

The mediating direction and the a-path recover cleanly; the b-path / indirect
effect is mildly attenuated at finite `T_i` (the per-subject covariance is
estimated from `T_i` rows).

## CAP-clustering — `capcov.clustering`

Clusters subjects whose covariance follows different covariate models, sharing a
projection, with membership informed by covariates.

```python
from capcov import clustering

# Y[i]: (T_i x p) array;  X: (n x q1) variance covariates (incl. intercept);
# W: (n x q2) membership covariates (incl. intercept)
fit = clustering.cap_pcl(Y, X, W, ncluster=2, stop_crt="nD", nD=1)

fit["gamma"]   # p x k shared projection loadings
fit["beta"]    # list of (q1 x K) per-cluster variance effects
fit["class"]   # n x k cluster assignment per direction

cf = clustering.cap_pcl_coef(Y, X, W, fit["gamma"][:, 0], ncluster=2)   # at a fixed gamma
cf["class"], cf["logLik"]
```

## Built-in examples

```python
from capcov import examples
examples.hdcap_example()       # -> {X, Y, truth}
examples.coc_example()         # -> {Y, X, W, truth}
examples.lcap_example()        # -> {Y, X, cov_names, truth}
examples.mediation_example()   # -> {X, M, Y, truth}
examples.clustering_example()  # -> {Y, X, W, truth}
```

Each `truth["gamma"]` (etc.) lets you check recovery, e.g. cosine similarity of a
fitted loading to its true direction.

## References

Get the citation(s) for any method programmatically:

```python
import capcov
capcov.references()        # all methods
capcov.references("coc")   # one method
capcov.REFERENCES          # the underlying dict
```
