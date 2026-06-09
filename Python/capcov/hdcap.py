"""HDCAP -- High-Dimensional Covariate Assisted Principal regression.

Faithful port of ``HDCAP/V6/CAP_HD.R`` (the de-sparsified version). With
``cov_shrinkage=False`` this is exactly **CAP** (Zhao et al. 2021); with
``cov_shrinkage=True`` it adds linear (Ledoit--Wolf) covariance shrinkage and is
auto-enabled when ``min_i T_i - 5 < p``.

Model:  ``log(gamma' Sigma_i gamma) = x_i' beta``.

Main entry point: :func:`cap_reg`.
"""
from __future__ import annotations

import numpy as np

from ._core import as_matrix_list, sigma_cube, scores, accum, gamma_solve, diag_level

__all__ = ["cap_reg", "cap_beta", "cap_beta_boot"]


# ---------------------------------------------------------------- weights -----
def _shrink_weights(score, Xbeta, Tvec, gg):
    """cap_beta-style linear-shrinkage weights (uses /Tvec and /(psi2+phi2))."""
    mu = np.mean(np.exp(Xbeta)) / gg
    delta2 = np.mean((score - mu * gg) ** 2)
    psi2 = min(np.mean((score - np.exp(Xbeta)) ** 2 / Tvec), delta2)
    phi2 = delta2 - psi2
    rho1 = psi2 * mu / (psi2 + phi2)
    rho2 = phi2 / (psi2 + phi2)
    return dict(rho1=rho1, rho2=rho2, mu=mu, phi2=phi2, psi2=psi2, delta2=delta2)


def _obj_weights(score, Xbeta, gg):
    """obj.func-style weights (no /Tvec, divides by delta2)."""
    mu = np.mean(np.exp(Xbeta)) / gg
    delta2 = np.mean((score - mu * gg) ** 2)
    psi2 = min(np.mean((score - np.exp(Xbeta)) ** 2), delta2)
    phi2 = delta2 - psi2
    return psi2 * mu / delta2, phi2 / delta2


def _obj_func(Sigma_raw, Tvec, X, gamma, beta):
    """CAP/HDCAP objective (always uses the raw sample covariance Sigma_raw)."""
    gamma = np.asarray(gamma, float).ravel()
    Xbeta = X @ beta
    sc = scores(Sigma_raw, gamma)
    gg = gamma @ gamma
    rho1, rho2 = _obj_weights(sc, Xbeta, gg)
    sk = rho1 * gg + rho2 * sc                      # gamma'(rho1 I + rho2 Sigma_i)gamma
    return float(np.sum(Tvec * (Xbeta + np.exp(-Xbeta) * sk)))


# ------------------------------------------------------- beta given gamma -----
def _cap_beta_core(Sigma, Tvec, X, gamma, cov_shrinkage, max_itr=1000, tol=1e-4):
    gamma = np.asarray(gamma, float).ravel()
    q = X.shape[1]
    sc = scores(Sigma, gamma)                       # raw projected variance
    gg = gamma @ gamma
    beta = np.zeros(q)
    rho1, rho2 = 0.0, 1.0
    if cov_shrinkage:
        w = _shrink_weights(sc, X @ beta, Tvec, gg)
        rho1, rho2 = w["rho1"], w["rho2"]
    s, diff = 0, np.inf
    while s <= max_itr and diff > tol:
        s += 1
        sk = (rho1 * gg + rho2 * sc) if cov_shrinkage else sc
        e = np.exp(-(X @ beta))
        wts = Tvec * sk * e
        Q1 = X.T @ (X * wts[:, None])
        Q2 = X.T @ (Tvec * (1.0 - sk * e))
        beta_new = beta - np.linalg.pinv(Q1) @ Q2
        if cov_shrinkage:
            w = _shrink_weights(sc, X @ beta_new, Tvec, gg)
            rho1, rho2 = w["rho1"], w["rho2"]
        diff = np.max(np.abs(beta_new - beta))
        beta = beta_new
    out = dict(beta=beta, convergence=(s < max_itr), score=sc)
    if cov_shrinkage:
        out["score_shrinkage"] = rho1 * gg + rho2 * sc
        out["shrinkage"] = _shrink_weights(sc, X @ beta, Tvec, gg)
    return out


def cap_beta(Y, X, gamma, cov_shrinkage=True, max_itr=1000, tol=1e-4):
    """Estimate ``beta`` for a fixed projection ``gamma`` (covariates ``X``)."""
    Y = as_matrix_list(Y)
    X = np.asarray(X, float)
    p = Y[0].shape[1]
    Tvec = np.array([Yi.shape[0] for Yi in Y], float)
    if Tvec.min() - 5 < p:
        cov_shrinkage = True
    Sigma, Tvec = sigma_cube(Y, cov_shrinkage=cov_shrinkage)
    return _cap_beta_core(Sigma, Tvec, X, gamma, cov_shrinkage, max_itr, tol)


# ------------------------------------------------- one direction (fixed pt) ---
def _cap_d1(Sigma, Tvec, X, H, gamma0, cov_shrinkage, max_itr=1000, tol=1e-4):
    """Alternating (beta, gamma) fixed point; returns the converged (unscaled) gamma."""
    n = Sigma.shape[2]
    p = Sigma.shape[0]
    gamma = np.asarray(gamma0, float).ravel()
    beta = np.zeros(X.shape[1])
    sc0 = scores(Sigma, gamma)                      # fixed initial-gamma score (R uses this in rho)
    gg = gamma @ gamma
    rho1, rho2 = 0.0, 1.0
    if cov_shrinkage:
        w = _shrink_weights(sc0, X @ beta, Tvec, gg)
        rho1, rho2 = w["rho1"], w["rho2"]
    s, diff = 0, np.inf
    while s <= max_itr and diff > tol:
        s += 1
        cur_gg = gamma @ gamma
        cur_sc = scores(Sigma, gamma)
        sk = (rho1 * cur_gg + rho2 * cur_sc) if cov_shrinkage else cur_sc
        e = np.exp(-(X @ beta))
        Q1 = X.T @ (X * (Tvec * sk * e)[:, None])
        Q2 = X.T @ (Tvec * (1.0 - sk * e))
        beta_new = beta - np.linalg.pinv(Q1) @ Q2
        if cov_shrinkage:
            # R uses the fixed initial-gamma score sc0 with the current gg
            w = _shrink_weights(sc0, X @ beta_new, Tvec, cur_gg)
            rho1, rho2 = w["rho1"], w["rho2"]
        e2 = np.exp(-(X @ beta_new))
        ww = Tvec * e2
        S1 = accum(Sigma, ww)
        if cov_shrinkage:
            S1 = rho1 * np.sum(ww) * np.eye(p) + rho2 * S1
        gamma_new = gamma_solve(S1, H)
        diff = max(np.max(np.abs(gamma_new - gamma)), np.max(np.abs(beta_new - beta)))
        gamma, beta = gamma_new, beta_new
    return gamma


def _default_gamma0_mat(p, seed=500, ncol=None):
    rng = np.random.RandomState(seed)
    G = rng.normal(size=((p + 1 + 5), p)).T          # p x (p+6)
    G = G / np.sqrt((G ** 2).sum(axis=0, keepdims=True))
    return G


def _cap_d1_opt(Y, X, cov_shrinkage, ninitial=None, seed=500, max_itr=1000, tol=1e-4):
    Y = as_matrix_list(Y)
    X = np.asarray(X, float)
    p = Y[0].shape[1]
    Sigma, Tvec = sigma_cube(Y, cov_shrinkage=cov_shrinkage)
    H = accum(Sigma, Tvec) / Tvec.sum()
    Sigma_raw, _ = sigma_cube(Y, cov_shrinkage=False)  # objective uses raw cov
    G = _default_gamma0_mat(p, seed=seed)
    k = G.shape[1] if ninitial is None else min(ninitial, G.shape[1])
    rng = np.random.RandomState(seed)
    cols = np.sort(rng.choice(G.shape[1], size=min(k, G.shape[1]), replace=False))
    G = G[:, cols]
    best, best_obj = None, np.inf
    for j in range(G.shape[1]):
        try:
            gu = _cap_d1(Sigma, Tvec, X, H, G[:, j], cov_shrinkage, max_itr, tol)
        except np.linalg.LinAlgError:
            continue
        bt = _cap_beta_core(Sigma, Tvec, X, gu, cov_shrinkage, max_itr, tol)
        obj = _obj_func(Sigma_raw, Tvec, X, gu, bt["beta"])
        if obj < best_obj:
            best_obj, best = obj, gu
    if best is None:
        raise RuntimeError("all initializations failed")
    # finalize: normalize + sign + refit beta on the normalized gamma
    gamma = best / np.sqrt(np.sum(best ** 2))
    if gamma[np.argmax(np.abs(gamma))] < 0:
        gamma = -gamma
    fin = _cap_beta_core(Sigma, Tvec, X, gamma, cov_shrinkage, max_itr, tol)
    out = dict(gamma=gamma, beta=fin["beta"], convergence=fin["convergence"],
               score=fin["score"])
    if cov_shrinkage:
        out["score_shrinkage"] = fin["score_shrinkage"]
        out["shrinkage"] = fin["shrinkage"]
    return out, H


def _cap_dk(Y, X, Phi0, cov_shrinkage, ninitial=None, seed=500, max_itr=1000, tol=1e-4):
    """Estimate the next direction after projecting out columns of ``Phi0``."""
    Y = as_matrix_list(Y)
    n = len(Y)
    p = Y[0].shape[1]
    p0 = Phi0.shape[1]
    Sigma, Tvec = sigma_cube(Y, cov_shrinkage=cov_shrinkage)
    beta0 = np.array([
        _cap_beta_core(Sigma, Tvec, X, Phi0[:, j], cov_shrinkage, max_itr, tol)["beta"][0]
        for j in range(p0)
    ])
    Ytmp = []
    for i, Yi in enumerate(Y):
        Y2 = Yi - Yi @ (Phi0 @ Phi0.T)
        if cov_shrinkage:
            Ytmp.append(Y2)                         # no add-back under shrinkage
        else:
            U, d, Vt = np.linalg.svd(Y2, full_matrices=False)
            ntmp = Y2.shape[0]
            dnew = np.concatenate([d[: (p - p0)], np.sqrt(np.exp(beta0) * ntmp)])
            Ytmp.append((U[:, : len(dnew)] * dnew) @ Vt[: len(dnew), :])
    res, _ = _cap_d1_opt(Ytmp, X, cov_shrinkage, ninitial, seed, max_itr, tol)
    res["orthogonal"] = Phi0.T @ res["gamma"]
    return res


# ------------------------------------------------------------- main entry -----
def cap_reg(Y, X, stop_crt="nD", nD=None, DfD_thred=5.0, cov_shrinkage=False,
            ninitial=None, seed=500, max_itr=1000, tol=1e-4, verbose=False):
    """Estimate CAP / HDCAP directions.

    Parameters
    ----------
    Y : list of ``T_i x p`` arrays (one per subject).
    X : ``n x q`` covariate matrix (include an intercept column).
    stop_crt : ``"nD"`` (fixed number) or ``"DfD"`` (deviation-from-diagonality).
    nD : number of directions when ``stop_crt="nD"``.
    DfD_thred : keep adding directions while DfD < this (when ``stop_crt="DfD"``).
    cov_shrinkage : linear covariance shrinkage; auto-enabled when min T_i - 5 < p.
                    ``False`` reproduces plain **CAP**.

    Returns a dict with ``gamma`` (p x nD), ``beta`` (q x nD), ``score``,
    ``orthogonality``, ``DfD`` (avg, sub), ``nD``, and (if shrinkage) ``shrinkage``.
    """
    Y = as_matrix_list(Y)
    X = np.asarray(X, float)
    n = len(Y)
    p = Y[0].shape[1]
    Tvec = np.array([Yi.shape[0] for Yi in Y], float)
    if Tvec.min() - 5 < p:
        cov_shrinkage = True
    if stop_crt == "nD" and nD is None:
        stop_crt = "DfD"

    re1, _ = _cap_d1_opt(Y, X, cov_shrinkage, ninitial, seed, max_itr, tol)
    Phi = re1["gamma"][:, None]
    beta = re1["beta"][:, None]
    score = re1["score"][:, None]
    sk_cols = [re1.get("shrinkage")]
    if verbose:
        print("Component 1")

    def _grow(stop_check):
        nonlocal Phi, beta, score
        while True:
            if not stop_check():
                break
            try:
                rk = _cap_dk(Y, X, Phi, cov_shrinkage, ninitial, seed, max_itr, tol)
            except Exception:
                break
            cand = np.column_stack([Phi, rk["gamma"]])
            if stop_crt == "DfD":
                avg, _ = diag_level(Y, cand)
                if avg[cand.shape[1] - 1] >= DfD_thred:
                    break
            Phi = cand
            beta = np.column_stack([beta, rk["beta"]])
            score = np.column_stack([score, rk["score"]])
            sk_cols.append(rk.get("shrinkage"))
            if verbose:
                print("Component", Phi.shape[1])
            if stop_crt == "nD" and Phi.shape[1] >= nD:
                break

    if stop_crt == "nD":
        if nD and nD > 1:
            _grow(lambda: Phi.shape[1] < nD)
    else:
        _grow(lambda: True)

    if Phi.shape[1] > 1:
        avg, sub = diag_level(Y, Phi)
        DfD = dict(avg_level=avg, sub_level=sub)
    else:
        DfD = dict(avg_level=np.array([1.0]), sub_level=np.ones((n, 1)))

    out = dict(gamma=Phi, beta=beta, score=score, nD=Phi.shape[1],
               orthogonality=Phi.T @ Phi, DfD=DfD, cov_shrinkage=cov_shrinkage)
    if cov_shrinkage and all(s is not None for s in sk_cols):
        out["shrinkage"] = sk_cols
    return out


# --------------------------------------------------------- bootstrap beta -----
def cap_beta_boot(Y, X, gamma, cov_shrinkage=True, sims=1000, conf_level=0.95,
                  seed=100, max_itr=1000, tol=1e-4):
    """Subject-level bootstrap inference for ``beta`` at a fixed ``gamma``."""
    Y = as_matrix_list(Y)
    X = np.asarray(X, float)
    n, q = len(Y), X.shape[1]
    boot = np.full((sims, q), np.nan)
    for b in range(sims):
        rng = np.random.RandomState(seed + b)
        idx = rng.randint(0, n, size=n)
        Yb = [Y[i] for i in idx]
        Xb = X[idx, :]
        try:
            res = cap_beta(Yb, Xb, gamma, cov_shrinkage=cov_shrinkage,
                           max_itr=max_itr, tol=tol)
            if res["convergence"]:
                boot[b, :] = res["beta"]
        except Exception:
            pass
    est = np.nanmean(boot, axis=0)
    se = np.nanstd(boot, axis=0, ddof=1)
    stat = est / se
    from scipy.stats import norm
    pval = 2 * (1 - norm.cdf(np.abs(stat)))
    z = norm.ppf(1 - (1 - conf_level) / 2)
    return dict(estimate=est, se=se, statistic=stat, pvalue=pval,
                lower=est - z * se, upper=est + z * se, boot=boot)
