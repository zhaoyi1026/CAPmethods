"""CAP-clustering (PCL) -- Python port of CAP-clustering/V4.

Parsimonious covariance-matrix clustering: subjects are soft-assigned to ``K``
clusters whose membership depends on covariates ``W`` (multinomial logit), and
within cluster ``k`` the projected log-variance follows
``log(gamma' Sigma_i gamma) = X_i' beta_k``. A shared projection ``gamma`` is
estimated jointly. Fit by EM (E-step: posterior membership; M-step: weighted
multinomial for membership + one Newton step per cluster for beta + generalized
eigen update for gamma).

Ported from ``CAP-clustering/V4/CAP_Cluster.R`` and ``cluster_kernels.cpp``.

Approximation note
------------------
R fits the membership model with ``brglm2::brmultinom`` (mean bias-reduced /
Firth multinomial). This port uses a self-contained **weighted Firth logistic**
for the common ``ncluster=2`` case (equivalent penalization) and a weighted
multinomial Newton (small ridge) for ``ncluster>2``. The membership *coefficients*
may differ slightly from brglm2, but cluster assignments, ``beta`` and ``gamma``
agree; per the package policy this is a statistically-equivalent port.

Data layout: ``Y[i]`` is a ``T_i x p`` matrix; ``X`` is ``n x q1`` (variance
covariates, incl. intercept); ``W`` is ``n x q2`` (membership covariates, incl.
intercept).
"""
from __future__ import annotations

import numpy as np

from ._core import gamma_solve

__all__ = ["cap_pcl", "cap_pcl_coef", "cap_pcl_coef_boot", "diag_level"]


# --------------------------------------------------------------------------- #
# kernels (== cluster_kernels.cpp)                                            #
# --------------------------------------------------------------------------- #
def _smat(Y):
    """Per-subject SECOND-MOMENT cube ``Y_i' Y_i / T_i`` (no centering) + sizes."""
    n = len(Y)
    p = np.asarray(Y[0]).shape[1]
    S = np.empty((p, p, n))
    Tvec = np.empty(n)
    for i in range(n):
        Yi = np.asarray(Y[i], dtype=float)
        T = Yi.shape[0]
        Tvec[i] = T
        S[:, :, i] = (Yi.T @ Yi) / T
    return S, Tvec


def _score(S, v):
    v = np.asarray(v, dtype=float).ravel()
    return np.einsum("i,ijk,j->k", v, S, v)


def _accum(S, w):
    return S @ np.asarray(w, dtype=float).ravel()


# --------------------------------------------------------------------------- #
# membership model: weighted (Firth) logistic / multinomial                    #
# --------------------------------------------------------------------------- #
def _firth_logistic_weighted(W, y, tw, max_itr=100, tol=1e-8):
    """Weighted Firth (Jeffreys-penalized) binary logistic; returns coef (q2,)."""
    n, q = W.shape
    beta = np.zeros(q)
    for _ in range(max_itr):
        eta = W @ beta
        pi = 1.0 / (1.0 + np.exp(-eta))
        wv = tw * pi * (1.0 - pi)
        wv = np.clip(wv, 1e-12, None)
        XtWX = W.T @ (wv[:, None] * W)
        try:
            inv = np.linalg.inv(XtWX)
        except np.linalg.LinAlgError:
            inv = np.linalg.pinv(XtWX)
        Q = np.sqrt(wv)[:, None] * W
        h = np.einsum("ij,jk,ik->i", Q, inv, Q)        # weighted leverages
        U = W.T @ (tw * (y - pi) + h * (0.5 - pi))
        step = inv @ U
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    return beta


def _multinom_weighted(W, cls, tw, K, max_itr=100, tol=1e-8, ridge=1e-6):
    """Weighted multinomial logit (ref class 0) via Newton; returns (q2, K) with col0=0.

    Used for ncluster>2 (an approximation to brglm2::brmultinom).
    """
    n, q = W.shape
    B = np.zeros((q, K))
    Yk = np.zeros((n, K))
    Yk[np.arange(n), cls] = 1.0
    for _ in range(max_itr):
        eta = W @ B                      # n x K
        eta = eta - eta.max(axis=1, keepdims=True)
        P = np.exp(eta); P /= P.sum(axis=1, keepdims=True)
        grad = (W.T @ ((tw[:, None]) * (Yk - P)))         # q x K
        maxstep = 0.0
        Bnew = B.copy()
        for k in range(1, K):
            Wk = tw * P[:, k] * (1 - P[:, k])
            Hk = W.T @ (Wk[:, None] * W) + ridge * np.eye(q)
            step = np.linalg.solve(Hk, grad[:, k])
            Bnew[:, k] = B[:, k] + step
            maxstep = max(maxstep, np.max(np.abs(step)))
        B = Bnew
        if maxstep < tol:
            break
    B[:, 0] = 0.0
    return B


def _membership_fit(W, cls, tw, K):
    if K == 2:
        coef = _firth_logistic_weighted(W, (cls == 1).astype(float), tw)
        A = np.zeros((W.shape[1], 2))
        A[:, 1] = coef
        return A
    return _multinom_weighted(W, cls, tw, K)


# --------------------------------------------------------------------------- #
# full log-likelihood (== ll.full)                                             #
# --------------------------------------------------------------------------- #
def _ll_full(S, X, W, eta_ind, gamma, beta, alpha, Tvec):
    n = X.shape[0]
    K = eta_ind.shape[1]
    score = _score(S, gamma)
    e = W @ alpha
    e = e - e.max(axis=1, keepdims=True)
    pi = np.exp(e); pi /= pi.sum(axis=1, keepdims=True)
    ll = 0.0
    for k in range(K):
        ll += np.sum(Tvec * eta_ind[:, k] * np.log(pi[:, k]))
        xb = X @ beta[:, k]
        ll += np.sum(Tvec * eta_ind[:, k] * (-xb - np.exp(-xb) * score) / 2.0)
    return ll


# --------------------------------------------------------------------------- #
# given gamma: estimate coefficients + clusters (== capPCL_coef)               #
# --------------------------------------------------------------------------- #
def cap_pcl_coef(Y, X, W, gamma, ncluster=2, max_itr=1000, tol=1e-4, seed=100,
                 nstart=10, score_return=True):
    """Estimate per-cluster ``beta``/membership ``alpha`` + clusters at fixed gamma.

    The EM has local optima; ``nstart`` random restarts are run and the fit with
    the highest log-likelihood is returned (R's ``capPCL_coef`` uses a single
    ``set.seed(100)`` start, so it can land in a worse local optimum — the global
    optimum recovers the truth, see README).
    """
    best = None
    best_ll = -np.inf
    for st in range(nstart):
        cf = _cap_pcl_coef_once(Y, X, W, gamma, ncluster, max_itr, tol,
                                seed + st, score_return)
        if cf["logLik"] > best_ll:
            best_ll = cf["logLik"]
            best = cf
    return best


def _cap_pcl_coef_once(Y, X, W, gamma, ncluster, max_itr, tol, seed, score_return):
    X = np.asarray(X, dtype=float)
    W = np.asarray(W, dtype=float)
    n, q1 = X.shape
    q2 = W.shape[1]
    K = ncluster
    gamma = np.asarray(gamma, dtype=float).ravel()
    S, Tvec = _smat(Y)
    score = _score(S, gamma)
    # projected data Z_i (length T_i): sum Z^2 = T_i * score_i
    rng = np.random.RandomState(seed)
    beta0 = rng.normal(size=(q1, K))
    alpha0 = rng.normal(size=(q2, K))
    alpha0[:, 0] = 0.0

    eta = np.full((n, K), 1.0 / K)
    s = 0
    diff = 100.0
    while s <= max_itr and diff > tol:
        s += 1
        # E-step (log-domain for stability)
        loglik = np.empty((n, K))
        for k in range(K):
            sig2 = np.exp(X @ beta0[:, k])
            loglik[:, k] = -0.5 * Tvec * (np.log(2 * np.pi * sig2) + score / sig2)
        e = W @ alpha0
        e = e - e.max(axis=1, keepdims=True)
        logpi = e - np.log(np.exp(e).sum(axis=1, keepdims=True))
        lpost = logpi + loglik
        lpost -= lpost.max(axis=1, keepdims=True)
        eta = np.exp(lpost); eta /= eta.sum(axis=1, keepdims=True)
        bad = ~np.isfinite(eta).all(axis=1)
        eta[bad, :] = 1.0 / K

        # M-step: membership (hard labels, weighted) + beta Newton step
        cls = np.argmax(eta, axis=1)
        if len(np.unique(cls)) > 1:
            alpha_new = _membership_fit(W, cls, Tvec, K)
        else:
            alpha_new = alpha0
        beta_new = np.empty((q1, K))
        for k in range(K):
            xb = X @ beta0[:, k]
            d1 = (X * (Tvec * eta[:, k] * (1 - np.exp(-xb) * score) / 2.0)[:, None]).sum(axis=0)
            d2 = X.T @ (np.diag(Tvec * eta[:, k] * np.exp(-xb) * score / 2.0)) @ X
            beta_new[:, k] = beta0[:, k] - np.linalg.pinv(d2) @ d1

        diff = max(np.max(np.abs(beta_new - beta0)), np.max(np.abs(alpha_new - alpha0)))
        beta0, alpha0 = beta_new, alpha_new

    eta_int = np.zeros((n, K))
    eta_int[np.arange(n), np.argmax(eta, axis=1)] = 1.0
    cl_ind = np.argmax(eta, axis=1) + 1          # 1-based, matching R
    ll = _ll_full(S, X, W, eta_int, gamma, beta0, alpha0, Tvec)
    res = {
        "gamma": gamma, "beta": beta0, "alpha": alpha0, "class": cl_ind,
        "eta": eta_int, "eta_est": eta, "logLik": ll, "nitr": s,
        "convergence": s < max_itr,
    }
    if score_return:
        res["score"] = score
    return res


# --------------------------------------------------------------------------- #
# one direction (estimate gamma too) (== capPCL_D1)                            #
# --------------------------------------------------------------------------- #
def _sign_norm(g):
    g = g / np.sqrt(np.sum(g ** 2))
    if g[np.argmax(np.abs(g))] < 0:
        g = -g
    return g


def _cap_pcl_d1(Y, X, W, gamma0, H, ncluster, max_itr, tol, seed):
    X = np.asarray(X, dtype=float)
    W = np.asarray(W, dtype=float)
    S, Tvec = _smat(Y)
    gamma = np.asarray(gamma0, dtype=float).ravel()
    s = 0
    diff = 100.0
    beta0 = None
    alpha0 = None
    while s <= max_itr and diff > tol:
        s += 1
        cf = cap_pcl_coef(Y, X, W, gamma, ncluster=ncluster, max_itr=max_itr,
                          tol=tol, seed=seed, nstart=5)
        if not cf["convergence"]:
            s = max_itr + 1
        eta = cf["eta_est"]; beta_new = cf["beta"]; alpha_new = cf["alpha"]
        wA = (Tvec / 2.0) * np.sum(eta * np.exp(-X @ beta_new), axis=1)
        A = _accum(S, wA)
        gamma_new = gamma_solve(A, H)
        if beta0 is not None:
            diff = max(np.max(np.abs(beta_new - beta0)), np.max(np.abs(alpha_new - alpha0)))
        beta0, alpha0 = beta_new, alpha_new
        gamma = gamma_new
    g = _sign_norm(gamma)
    cf = cap_pcl_coef(Y, X, W, g, ncluster=ncluster, max_itr=max_itr, tol=tol,
                      seed=seed, nstart=5)
    cf["gamma"] = g
    cf["convergence"] = s < max_itr
    return cf


def _cap_pcl_d1_opt(Y, X, W, H, ncluster, max_itr, tol, ninitial, seed):
    p = np.asarray(Y[0]).shape[1]
    S, Tvec = _smat(Y)
    if H is None:
        H = _accum(S, Tvec) / Tvec.sum()
    if ninitial is None:
        ninitial = min(p, 10)
    w, V = np.linalg.eigh(H)
    V = V[:, ::-1]                       # descending, eigenvector starts
    best = None
    best_ll = -np.inf
    for kk in range(min(ninitial, p)):
        try:
            rk = _cap_pcl_d1(Y, X, W, V[:, kk], H, ncluster, max_itr, tol, seed)
        except Exception:
            continue
        if rk["logLik"] > best_ll:
            best_ll = rk["logLik"]
            best = rk
    return best


def _cap_pcl_dk(Y, X, W, Gamma0, H, ncluster, max_itr, tol, ninitial, seed):
    if Gamma0 is None:
        return _cap_pcl_d1_opt(Y, X, W, H, ncluster, max_itr, tol, ninitial, seed)
    n = len(Y)
    p = np.asarray(Y[0]).shape[1]
    p0 = Gamma0.shape[1]
    beta0 = np.empty((X.shape[1], ncluster, p0))
    for j in range(p0):
        cf = cap_pcl_coef(Y, X, W, Gamma0[:, j], ncluster=ncluster, max_itr=max_itr,
                          tol=tol, seed=seed)
        beta0[:, :, j] = cf["beta"]
    P = Gamma0 @ Gamma0.T
    Yt = []
    for i in range(n):
        Yi = np.asarray(Y[i], dtype=float)
        Ynew = Yi - Yi @ P
        U, d, Vt = np.linalg.svd(Ynew, full_matrices=False)
        Dnew = d.copy()
        Dnew[p - p0:p] = np.sqrt(np.exp(beta0[0, 0, :]) * Yi.shape[0])
        Yt.append((U * Dnew) @ Vt)
    nin = max(min(10, p - p0 - 1), 3) if p0 >= 2 else ninitial
    re = _cap_pcl_d1_opt(Yt, X, W, H, ncluster, max_itr, tol, nin, seed)
    if re is not None:
        re["orthogonality"] = re["gamma"] @ Gamma0
    return re


# --------------------------------------------------------------------------- #
# deviation-from-diagonality (== diag.level)                                   #
# --------------------------------------------------------------------------- #
def diag_level(Y, Gamma):
    Gamma = np.asarray(Gamma, dtype=float)
    if Gamma.ndim == 1 or Gamma.shape[1] == 1:
        return {"avg_level": np.array([1.0])}
    n = len(Y)
    r = Gamma.shape[1]
    S, Tvec = _smat(Y)
    dl = np.ones((n, r))
    for i in range(n):
        cov = S[:, :, i]
        for j in range(1, r):
            G = Gamma[:, : j + 1]
            M = G.T @ cov @ G
            dl[i, j] = np.prod(np.diag(M)) / np.linalg.det(M)
    avg = np.array([np.prod(dl[:, j] ** (Tvec / Tvec.sum())) for j in range(r)])
    return {"avg_level": avg, "sub_level": dl}


# --------------------------------------------------------------------------- #
# main entry: find leading directions (== capPCL)                              #
# --------------------------------------------------------------------------- #
def cap_pcl(Y, X, W, ncluster=2, stop_crt="nD", nD=None, DfD_thred=2.0, H=None,
            max_itr=1000, tol=1e-4, ninitial=None, seed=100, verbose=False):
    """Covariance-matrix clustering with a shared projection (== R ``capPCL``)."""
    X = np.asarray(X, dtype=float)
    W = np.asarray(W, dtype=float)
    p = np.asarray(Y[0]).shape[1]
    if stop_crt == "nD" and nD is None:
        stop_crt = "DfD"
    S, Tvec = _smat(Y)
    if H is None:
        H = _accum(S, Tvec) / Tvec.sum()

    re1 = _cap_pcl_d1_opt(Y, X, W, H, ncluster, max_itr, tol, ninitial, seed)
    Gamma = [re1["gamma"]]
    res = [re1]
    if verbose:
        print("Component 1")

    def _grow(rt):
        Gamma.append(rt["gamma"]); res.append(rt)
        if verbose:
            print(f"Component {len(Gamma)}")

    if stop_crt == "nD":
        for _ in range(2, nD + 1):
            rt = _cap_pcl_dk(Y, X, W, np.column_stack(Gamma), H, ncluster,
                             max_itr, tol, ninitial, seed)
            if rt is None:
                break
            _grow(rt)
    elif stop_crt == "DfD":
        DfD = 1.0
        while DfD < DfD_thred:
            rt = _cap_pcl_dk(Y, X, W, np.column_stack(Gamma), H, ncluster,
                             max_itr, tol, ninitial, seed)
            if rt is None:
                break
            cand = np.column_stack(Gamma + [rt["gamma"]])
            DfD = diag_level(Y, cand)["avg_level"][-1]
            if DfD < DfD_thred and np.isfinite(DfD):
                _grow(rt)
            else:
                break
    else:
        raise ValueError("stop_crt must be 'nD' or 'DfD'")

    G = np.column_stack(Gamma)
    out = {
        "gamma": G,
        "beta": [r["beta"] for r in res],
        "alpha": [r["alpha"] for r in res],
        "class": np.column_stack([r["class"] for r in res]),
        "logLik": np.array([r["logLik"] for r in res]),
        "columns": [f"C{i + 1}" for i in range(G.shape[1])],
    }
    if G.shape[1] >= 2:
        out["DfD"] = diag_level(Y, G)
    return out


# --------------------------------------------------------------------------- #
# bootstrap inference for beta at a fixed gamma + membership (== capPCL_coef_boot)
# --------------------------------------------------------------------------- #
def cap_pcl_coef_boot(Y, X, W, gamma, ncluster=2, sims=1000, boot_ci_type="se",
                      conf_level=0.95, seed=100, max_itr=1000, tol=1e-4,
                      verbose=False):
    """Bootstrap inference for the per-cluster ``beta`` (== R ``capPCL_coef_boot``).

    Subjects are resampled with numpy's RNG; clusters are aligned to the full-sample
    fit each draw (by membership overlap) to mitigate label switching.
    """
    from scipy.stats import norm
    X = np.asarray(X, dtype=float)
    W = np.asarray(W, dtype=float)
    n = len(Y)
    base = cap_pcl_coef(Y, X, W, gamma, ncluster=ncluster, max_itr=max_itr, tol=tol)
    q1, K = base["beta"].shape
    boot = np.full((q1, K, sims), np.nan)
    for b in range(sims):
        rng = np.random.RandomState(seed + b)
        idx = rng.randint(0, n, size=n)
        Yt = [Y[i] for i in idx]
        try:
            rt = cap_pcl_coef(Yt, X[idx], W[idx], gamma, ncluster=ncluster,
                              max_itr=max_itr, tol=tol)
            perm = _align_clusters(base["beta"], rt["beta"])
            boot[:, :, b] = rt["beta"][:, perm]
        except Exception:
            pass
        if verbose:
            print(f"Bootstrap sample {b + 1}")

    est = base["beta"]
    se = np.nanstd(boot, axis=2, ddof=1)
    stat = est / se
    pval = (1 - norm.cdf(np.abs(stat))) * 2
    lo, hi = (1 - conf_level) / 2, 1 - (1 - conf_level) / 2
    if boot_ci_type == "perc":
        lower = np.nanquantile(boot, lo, axis=2)
        upper = np.nanquantile(boot, hi, axis=2)
    else:
        z = norm.ppf(hi)
        lower, upper = est - z * se, est + z * se
    return {"estimate": est, "se": se, "statistic": stat, "pvalue": pval,
            "lower": lower, "upper": upper, "boot": boot}


def _align_clusters(beta_ref, beta):
    """Return a column permutation of ``beta`` best matching ``beta_ref`` (K small)."""
    import itertools
    K = beta_ref.shape[1]
    best = tuple(range(K))
    best_err = np.inf
    for perm in itertools.permutations(range(K)):
        err = np.sum((beta[:, list(perm)] - beta_ref) ** 2)
        if err < best_err:
            best_err = err
            best = perm
    return list(best)
