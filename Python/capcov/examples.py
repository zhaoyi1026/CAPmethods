"""Built-in example-data generators, one per implemented method.

Each returns a dict with the inputs the matching wrapper expects plus a
``truth`` entry holding the data-generating parameters, so estimates can be
checked. These mirror the R ``CAPmethods`` ``*_example()`` generators (kept in
sync with that package; settings updated 2026-06-09 to the manuscript / Shiny-demo
simulations). Several use the published simulation settings:

- ``hdcap_example``     -- HD-shrinkage CAP: p=20, two covariate-driven directions.
- ``lcap_example``      -- longitudinal CAP: p=20, two within-subject covariates.
- ``coc_example``       -- case-1 sim of Zhao & Zhao (Biometrics 2025).
- ``mediation_example`` -- single-treatment covariance-mediator sim (Xu & Zhao,
  Biostatistics 2025), IE = alpha*beta = 0.8*0.7 = 0.56.
- ``clustering_example``-- p=50, two covariate-driven components (D2, D4).

Random data differ from the R run (different RNG), but the structure, the true
parameters, and the recovery behavior match. Only methods implemented in Python
have a generator here (HDCAP/CAP, CAP-CoC, LCAP, CAP-mediation, CAP-clustering);
MCAP and CAP-HDcov are still stubs.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "hdcap_example", "lcap_example", "coc_example",
    "mediation_example", "clustering_example",
]


def _root_sig(Phi, eigs):
    """Symmetric PSD square root of ``Phi diag(eigs) Phi'``."""
    eigs = np.asarray(eigs, dtype=float)
    S = Phi @ np.diag(eigs) @ Phi.T
    S = (S + S.T) / 2
    w, V = np.linalg.eigh(S)
    return V @ np.diag(np.sqrt(np.clip(w, 1e-8, None))) @ V.T


def _ortho_pos(M):
    """QR-orthonormalize and sign-fix each column to have positive max-|entry|."""
    Q, _ = np.linalg.qr(M)
    for j in range(Q.shape[1]):
        if Q[np.argmax(np.abs(Q[:, j])), j] < 0:
            Q[:, j] = -Q[:, j]
    return Q


def _recycled_basis(p, seed):
    """Orthonormal basis from ``runif(p)`` recycled into a p x p matrix (== R)."""
    v = np.random.RandomState(seed).uniform(size=p)
    return _ortho_pos(np.tile(v[:, None], (1, p)))


# --------------------------------------------------------------------------- #
# HDCAP / CAP -- HD-shrinkage sim, two covariate-driven directions             #
# --------------------------------------------------------------------------- #
def hdcap_example(n=100, Ti=100, seed=100):
    """Example for ``capcov.cap_reg`` / ``cap_beta`` (HDCAP / CAP).

    p=20 with a common eigenbasis and one binary covariate (group). Two directions
    (basis cols 2 and 4) satisfy the CAP model with group effects -1 and +1; every
    other direction has a covariate-free random log-variance. With covariance
    shrinkage on, ``nD=2`` recovers both. Returns ``X`` (n x 2), ``Y`` (list of
    T_i x p), and ``truth`` (gamma p x 2, beta 2 x 2).
    """
    p = 20
    beta_sd = 0.5
    Gamma = _recycled_basis(p, seed)
    beta1 = np.concatenate([np.linspace(5, 1, 10), np.linspace(0.5, -1, p - 10)])
    beta2 = np.concatenate([[0.0, -1.0, 0.0, 1.0, 0.0], np.zeros(p - 5)])
    beta_mat = np.vstack([beta1, beta2])               # 2 x p (Intercept, group)
    cov_dirs = np.where(beta2 != 0)[0]                 # -> [1, 3] (cols 2,4)

    rng = np.random.RandomState(seed)
    X = np.column_stack([np.ones(n), rng.binomial(1, 0.5, size=n)])
    Y = []
    for i in range(n):
        eigs = np.empty(p)
        for j in range(p):
            if beta2[j] == 0:
                eigs[j] = np.exp(rng.normal(beta1[j], beta_sd))
            else:
                eigs[j] = np.exp(X[i] @ beta_mat[:, j])
        Y.append(rng.normal(size=(Ti, p)) @ _root_sig(Gamma, eigs))
    return {"X": X, "Y": Y,
            "truth": {"gamma": Gamma[:, cov_dirs], "beta": beta_mat[:, cov_dirs]}}


# --------------------------------------------------------------------------- #
# LCAP -- longitudinal, two within-subject covariates, two directions          #
# --------------------------------------------------------------------------- #
def lcap_example(n=100, nV=5, Ti=100, seed=4):
    """Example for ``capcov.lcap.cap_reg`` / ``cap_beta``.

    p=20 with a common time-invariant eigenbasis, ~nV visits/subject, two
    within-subject covariates (x1, x2) and a subject random intercept. Two
    directions (basis cols 2 and 4) satisfy the CAP model; ``nD=2`` with shrinkage
    recovers both. Returns ``Y`` (list over subjects of lists of T_iv x p),
    ``X`` (list of nV_i x 3 designs), and ``truth`` (gamma p x 2, beta 3 x 2).
    """
    p = 20
    Gamma = _recycled_basis(p, 100)
    beta_mat = np.vstack([
        np.concatenate([np.linspace(3, -1, 5), np.linspace(-1.5, -3, p - 5)]),   # intercept
        np.concatenate([[0.0, -0.5, 0.0, 0.5, 0.0], np.zeros(p - 5)]),           # x1
        np.concatenate([[0.0, 0.5, 0.0, -0.25, 0.0], np.zeros(p - 5)]),          # x2
    ])
    cov_dirs = np.where(beta_mat[1] != 0)[0]           # -> [1, 3] (cols 2,4)

    rng = np.random.RandomState(seed)
    nVvec = np.maximum(2, np.round(rng.normal(nV, 1.0, size=n)).astype(int))
    beta0 = np.empty((n, p))
    for j in range(p):
        beta0[:, j] = rng.normal(beta_mat[0, j], 0.1, size=n)
    Y, X = [], []
    for i in range(n):
        nv = nVvec[i]
        Xi = np.column_stack([np.ones(nv), rng.normal(0, 0.5, nv), rng.normal(0, 0.5, nv)])
        vis = []
        for v in range(nv):
            Tiv = max(p + 5, int(round(rng.normal(Ti, 5))))
            delta = np.array([np.exp(Xi[v] @ np.concatenate([[beta0[i, j]], beta_mat[1:, j]]))
                              for j in range(p)])
            vis.append(rng.normal(size=(Tiv, p)) @ _root_sig(Gamma, delta))
        Y.append(vis); X.append(Xi)
    return {"Y": Y, "X": X, "cov_names": ["x1", "x2"],
            "truth": {"gamma": Gamma[:, cov_dirs], "beta": beta_mat[:, cov_dirs]}}


# --------------------------------------------------------------------------- #
# CAP-CoC -- case-1 sim of Zhao & Zhao (Biometrics 2025)                        #
# --------------------------------------------------------------------------- #
def coc_example(n=150, p=10, q=5, Tx=150, Ty=150, seed=2024):
    """Example for ``capcov.coc_reg`` (covariance-on-covariance).

    Case-1 simulation of Zhao & Zhao (Biometrics 2025): predictor directions
    {1, 3} drive outcome directions {2, 4} with ``alpha = (3, 2)`` and covariate
    effects ``beta``; the other directions carry covariate-free noise. Returns
    ``Y``, ``X`` (lists), ``W`` (n x 2), and ``truth`` (gamma, theta, alpha, beta).
    """
    if p < 3 or q < 4:
        raise ValueError("need p >= 3 and q >= 4")
    Gamma1 = _ortho_pos(np.random.RandomState(100).uniform(size=(p, p)))
    Gamma2 = _ortho_pos(np.random.RandomState(500).uniform(size=(q, q)))
    x_eig_m = np.exp(np.linspace(1, -2, p)); x_eig_sd = 0.5
    y_eig_m = np.exp(np.linspace(1, -2, q)); y_eig_sd = 0.5
    x_idx = [0, 2]; y_idx = [1, 3]              # 0-based ({1,3} and {2,4} in R)
    alpha = np.array([3.0, 2.0]); beta = np.array([[1.0, -1.0], [-1.0, 1.0]])

    rng = np.random.RandomState(seed)
    W = np.column_stack([np.ones(n), rng.binomial(1, 0.5, size=n)])
    L1 = np.empty((n, p))
    for j in range(p):
        L1[:, j] = np.exp(rng.normal(np.log(x_eig_m[j]), x_eig_sd, size=n))
    L2 = np.empty((n, q))
    for k in range(q):
        if k in y_idx:
            f = y_idx.index(k)
            L2[:, k] = np.exp(alpha[f] * np.log(L1[:, x_idx[f]]) + W @ beta[:, f])
        else:
            L2[:, k] = np.exp(rng.normal(np.log(y_eig_m[k]), y_eig_sd, size=n))
    X, Y = [], []
    for i in range(n):
        SX = Gamma1 @ np.diag(L1[i]) @ Gamma1.T
        SY = Gamma2 @ np.diag(L2[i]) @ Gamma2.T
        X.append(rng.multivariate_normal(np.zeros(p), SX, size=Tx))
        Y.append(rng.multivariate_normal(np.zeros(q), SY, size=Ty))
    return {"Y": Y, "X": X, "W": W,
            "truth": {"gamma": Gamma2[:, y_idx], "theta": Gamma1[:, x_idx],
                      "alpha": alpha, "beta": beta,
                      "recovering_settings": dict(Hy=np.eye(q), Hx=np.eye(p),
                                                  burn_in=200, nD=1)}}


# --------------------------------------------------------------------------- #
# CAP-mediation -- single-treatment covariance-mediator sim                     #
# --------------------------------------------------------------------------- #
def mediation_example(n=100, p=10, Ti=150, seed=2024):
    """Example for ``capcov.cap_med`` (covariance / graph mediator).

    Single-treatment simulation (Xu & Zhao, Biostatistics 2025): a binary exposure
    shifts the variance of a p-dim mediator along one latent direction ``theta``;
    the outcome depends on that log-variance (b-path ``beta=0.7``) plus a direct
    effect. Exposure -> mediator-variance is the a-path (``alpha=0.8``); the true
    indirect effect is ``IE = alpha*beta = 0.56``. Returns ``X`` (n x 1), ``M``
    (list of T_i x p), ``Y`` (n,), and ``truth`` (theta, alpha, beta, gamma, IE).
    """
    rng = np.random.RandomState(seed)
    alpha0, alpha_x, beta, gamma0, gamma_x = 0.2, 0.8, 0.7, 0.1, 0.4
    x = np.tile([0.0, 1.0], n // 2 + 1)[:n]
    X = x.reshape(-1, 1)
    Phi = _ortho_pos(rng.normal(size=(p, p)))
    base = np.concatenate([[np.nan], np.exp(np.linspace(np.log(1.2), np.log(0.3), p - 1))])
    M, Y = [], np.empty(n)
    for i in range(n):
        lv = alpha0 + alpha_x * x[i] + rng.normal(0, 0.3)
        ev = base.copy(); ev[0] = np.exp(lv)
        M.append(rng.normal(size=(Ti, p)) @ _root_sig(Phi, ev))
        Y[i] = gamma0 + gamma_x * x[i] + beta * lv + rng.normal(0, 0.3)
    return {"X": X, "M": M, "Y": Y,
            "truth": {"theta": Phi[:, 0], "alpha": alpha_x, "beta": beta,
                      "gamma": gamma_x, "IE": alpha_x * beta}}


# --------------------------------------------------------------------------- #
# CAP-clustering (PCL) -- p=50, two covariate-driven components (D2, D4)        #
# --------------------------------------------------------------------------- #
def clustering_example(n=100, p=50, Ti=100, seed=1):
    """Example for ``capcov.clustering.cap_pcl``.

    p=50, K=2 clusters, with two covariate-driven components -- D2 (well separated)
    and D4 (clusters differ only in covariate signs). Within a cluster the
    projected log-variance is a log-linear model in ``X``; membership depends on
    ``W``. Returns ``Y`` (list of T_i x p), ``X`` (n x 3), ``W`` (n x 2), and
    ``truth`` (gamma p x 2 [D2, D4], cluster n x 2, beta {D2, D4}).
    """
    rng = np.random.RandomState(seed)
    Pi = _ortho_pos(rng.normal(size=(p, p)))
    x1 = rng.binomial(1, 0.5, size=n); x2 = rng.normal(size=n)
    X = np.column_stack([np.ones(n), x1, x2])              # q1 = 3
    w1 = rng.binomial(1, 0.5, size=n)
    W = np.column_stack([np.ones(n), w1])                  # q2 = 2
    cD2 = 1 + rng.binomial(1, 1.0 / (1.0 + np.exp(-(W @ np.array([0.5, -1.0])))))
    cD4 = 1 + rng.binomial(1, 1.0 / (1.0 + np.exp(-(W @ np.array([-0.25, 0.5])))))
    bD2 = np.array([[1.0, -1.0], [1.0, -1.0], [-1.0, 1.0]])   # q1 x K (well separated)
    bD4 = np.array([[0.5, 0.5], [0.5, -0.5], [-0.5, 0.5]])    # same baseline
    logmean = np.linspace(3, -1, p)
    Y = []
    for i in range(n):
        lam = np.exp(rng.normal(logmean, 0.2))
        lam[1] = np.exp(X[i] @ bD2[:, cD2[i] - 1])
        lam[3] = np.exp(X[i] @ bD4[:, cD4[i] - 1])
        Y.append(rng.normal(size=(Ti, p)) @ _root_sig(Pi, lam))
    return {"Y": Y, "X": X, "W": W,
            "truth": {"gamma": np.column_stack([Pi[:, 1], Pi[:, 3]]),
                      "cluster": np.column_stack([cD2, cD4]),
                      "beta": {"D2": bD2, "D4": bD4}}}
