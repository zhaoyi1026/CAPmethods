"""HDCAP -- High-Dimensional Covariate Assisted Principal regression.

Port of the R ``capReg`` of HDCAP V9 (CAPmethods 0.2.0). With
``cov_shrinkage=False`` this is **CAP** (Zhao et al. 2021); with
``cov_shrinkage=True`` it adds linear (Ledoit--Wolf) covariance shrinkage and is
auto-enabled when ``min_i T_i - 5 < p``.

Model:  ``log(gamma' Sigma_i gamma) = x_i' beta``.

Main entry point: :func:`cap_reg`.

As HDCAP V9 (capcov 0.2.0):

* The objective that compares the initial values of a direction is the function
  the fit minimizes, ``sum_i T_i (x_i'beta + s_i exp(-x_i'beta))`` with
  ``gamma'H gamma = 1``: ``s_i = gamma'S_i gamma`` without shrinkage (the
  negative log-likelihood); with shrinkage ``s_i = rho1 gamma'gamma +
  rho2 gamma'Sigma_i gamma``, the shrunk score of the fit (``Sigma_i`` the
  per-subject ``cov_ls``, ``rho1, rho2`` the fit's weights at beta).
* Several initial values of beta, argument ``beta0``: a vector of length q or a
  ``q x m`` array (one per column); ``None`` = two, ``0`` and
  :func:`cap_beta0_v5` (R's ``set.seed(100); c(0, rnorm(q-1, 1, 1))``, the same
  numbers as in R).  Each initial direction runs from each initial value of
  beta; every beta fit keeps the converged estimate with the lowest objective.
"""
from __future__ import annotations

import numpy as np

from ._core import (as_matrix_list, sigma_cube, scores, accum, gamma_solve,
                    diag_level, r_rnorm)

__all__ = ["cap_reg", "cap_beta", "cap_beta_boot", "obj_func", "cap_beta0_v5"]


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


# ------------------------------------------------------- initial values ------
def cap_beta0_v5(q):
    """HDCAP's V5-V8 initial value of beta: R's ``set.seed(100); c(0, rnorm(q-1, 1, 1))``.

    The same numbers as in R (via :func:`capcov._core.r_rnorm`); numpy's random
    numbers are not changed.
    """
    return np.r_[0.0, r_rnorm(int(q) - 1, 100, 1.0, 1.0)]


def _beta0_mat(q, beta0=None):
    """Initial values of beta as a ``q x m`` array, one per column.

    ``None``: two, ``0`` and :func:`cap_beta0_v5`; a vector of length q: that
    one; a ``q x m`` array: its columns.
    """
    if beta0 is None:
        return np.column_stack([np.zeros(q), cap_beta0_v5(q)])
    B0 = np.array(beta0, dtype=float)
    if B0.ndim == 1:
        B0 = B0[:, None]
    if B0.ndim != 2 or B0.shape[0] != q or B0.shape[1] < 1 or not np.all(np.isfinite(B0)):
        raise ValueError(f"beta0 must be a finite vector of length ncol(X) = {q}, "
                         f"or an array with {q} rows (one initial value per column).")
    return B0


# -------------------------------------------------------------- objective ----
def _cap_obj(sc, gg, X, Tvec, beta, cov_shrinkage):
    """The function the fit minimizes: sum_i T_i (x_i'beta + s_i exp(-x_i'beta)).

    ``sc``: gamma'Sigma_i gamma with the covariances of the fit; ``s = sc``
    without shrinkage, ``rho1 gg + rho2 sc`` (the fit's weights at beta) with.
    """
    Xbeta = X @ beta
    s = sc
    if cov_shrinkage:
        w = _shrink_weights(sc, Xbeta, Tvec, gg)
        s = w["rho1"] * gg + w["rho2"] * sc
    return float(np.sum(Tvec * (Xbeta + np.exp(-Xbeta) * s)))


def obj_func(Y, X, gamma, beta, cov_shrinkage=False):
    """Objective of (gamma, beta) (== R ``obj.func`` of HDCAP V9).

    ``sum_i T_i (x_i'beta + s_i exp(-x_i'beta))`` with the covariances of the
    fit; ``cov_shrinkage`` is forced, as in the fit, when ``min_i T_i - 5 < p``.
    """
    Y = as_matrix_list(Y)
    X = np.asarray(X, float)
    p = Y[0].shape[1]
    if min(Yi.shape[0] for Yi in Y) - 5 < p:
        cov_shrinkage = True
    Sigma, Tvec = sigma_cube(Y, cov_shrinkage=cov_shrinkage)
    gamma = np.asarray(gamma, float).ravel()
    return _cap_obj(scores(Sigma, gamma), float(gamma @ gamma), X, Tvec,
                    np.asarray(beta, float).ravel(), cov_shrinkage)


# ------------------------------------------------------- beta given gamma -----
def _cap_beta_core(sc, gg, X, Tvec, cov_shrinkage, beta, max_itr=1000, tol=1e-4):
    """Newton iteration for beta from the initial value ``beta``; returns (beta, iter)."""
    beta = np.array(beta, dtype=float)
    rho1, rho2 = 0.0, 1.0
    if cov_shrinkage:
        w = _shrink_weights(sc, X @ beta, Tvec, gg)
        rho1, rho2 = w["rho1"], w["rho2"]
    s, diff = 0, np.inf
    while s <= max_itr and diff > tol:
        s += 1
        sk = (rho1 * gg + rho2 * sc) if cov_shrinkage else sc
        e = np.exp(-(X @ beta))
        Q1 = X.T @ (X * (Tvec * sk * e)[:, None])
        Q2 = X.T @ (Tvec * (1.0 - sk * e))
        beta_new = beta - np.linalg.pinv(Q1) @ Q2
        if cov_shrinkage:
            w = _shrink_weights(sc, X @ beta_new, Tvec, gg)
            rho1, rho2 = w["rho1"], w["rho2"]
        diff = np.max(np.abs(beta_new - beta))
        beta = beta_new
    return beta, s


def _cap_beta_fit(sc, gg, X, Tvec, cov_shrinkage, B0, max_itr=1000, tol=1e-4):
    """beta from every initial value (columns of ``B0``).

    Kept: the converged estimate with the lowest objective (a non-converged one
    only when none converged).  A start whose iteration fails or diverges is
    skipped; an error when none gives a finite estimate.
    """
    best = None
    for j in range(B0.shape[1]):
        try:
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                b, it = _cap_beta_core(sc, gg, X, Tvec, cov_shrinkage, B0[:, j], max_itr, tol)
                if not np.all(np.isfinite(b)):
                    continue
                ob = _cap_obj(sc, gg, X, Tvec, b, cov_shrinkage)
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            continue
        if not np.isfinite(ob):
            ob = np.inf
        conv = it < max_itr
        if (best is None or (conv and not best["conv"])
                or (conv == best["conv"] and ob < best["obj"])):
            best = dict(beta=b, conv=conv, obj=ob, start=j)
    if best is None:
        raise RuntimeError("cap_beta: beta could not be estimated from any initial value of beta.")
    return best


def _beta_result(sc, gg, X, Tvec, cov_shrinkage, fit):
    out = dict(beta=fit["beta"], convergence=fit["conv"], score=sc)
    if cov_shrinkage:
        w = _shrink_weights(sc, X @ fit["beta"], Tvec, gg)
        out["score_shrinkage"] = w["rho1"] * gg + w["rho2"] * sc
        out["shrinkage"] = w
    return out


def cap_beta(Y, X, gamma, cov_shrinkage=True, max_itr=1000, tol=1e-4, beta0=None):
    """Estimate ``beta`` for a fixed projection ``gamma`` (covariates ``X``).

    ``beta0``: initial value(s) of beta (see :func:`cap_reg`); the converged
    estimate with the lowest objective is kept.
    """
    Y = as_matrix_list(Y)
    X = np.asarray(X, float)
    p = Y[0].shape[1]
    Tvec = np.array([Yi.shape[0] for Yi in Y], float)
    if Tvec.min() - 5 < p:
        cov_shrinkage = True
    Sigma, Tvec = sigma_cube(Y, cov_shrinkage=cov_shrinkage)
    gamma = np.asarray(gamma, float).ravel()
    sc, gg = scores(Sigma, gamma), float(gamma @ gamma)
    fit = _cap_beta_fit(sc, gg, X, Tvec, cov_shrinkage, _beta0_mat(X.shape[1], beta0), max_itr, tol)
    return _beta_result(sc, gg, X, Tvec, cov_shrinkage, fit)


# ------------------------------------------------- one direction (fixed pt) ---
def _cap_d1(Sigma, Tvec, X, H, gamma0, cov_shrinkage, beta0, max_itr=1000, tol=1e-4):
    """Alternating (beta, gamma) fixed point from (gamma0, beta0); returns the
    converged (unscaled) gamma."""
    p = Sigma.shape[0]
    gamma = np.asarray(gamma0, float).ravel()
    beta = np.array(beta0, dtype=float)
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


def _sign_unit(g):
    g = np.asarray(g, float) / np.sqrt(np.sum(np.asarray(g, float) ** 2))
    return -g if g[np.argmax(np.abs(g))] < 0 else g


def _cap_dir_obj(gamma, Sigma, H, X, Tvec, cov_shrinkage, B0, max_itr, tol):
    """Objective of a direction, as the initial values are compared: gamma scaled
    to gamma'H gamma = 1, beta from every initial value of beta."""
    gu = gamma / np.sqrt(gamma @ H @ gamma)
    sc = scores(Sigma, gu)
    return _cap_beta_fit(sc, float(gu @ gu), X, Tvec, cov_shrinkage, B0, max_itr, tol)["obj"]


def _cap_d1_opt(Y, X, cov_shrinkage, ninitial=None, seed=500, max_itr=1000, tol=1e-4,
                beta0=None, gamma0_mat=None):
    Y = as_matrix_list(Y)
    X = np.asarray(X, float)
    p = Y[0].shape[1]
    Sigma, Tvec = sigma_cube(Y, cov_shrinkage=cov_shrinkage)
    H = accum(Sigma, Tvec) / Tvec.sum()
    G = _default_gamma0_mat(p, seed=seed) if gamma0_mat is None else np.array(gamma0_mat, float)
    if G.ndim == 1:
        G = G[:, None]
    k = min(G.shape[1], 10) if ninitial is None else min(int(ninitial), G.shape[1])   # R: min(ncol, 10)
    if k < G.shape[1]:
        rng = np.random.RandomState(seed)
        G = G[:, np.sort(rng.choice(G.shape[1], size=k, replace=False))]
    # every initial value of gamma from every initial value of beta
    B0 = _beta0_mat(X.shape[1], beta0)
    cands, objs = [], []
    for j in range(G.shape[1]):
        for jj in range(B0.shape[1]):
            try:
                with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                    gu = _cap_d1(Sigma, Tvec, X, H, G[:, j], cov_shrinkage, B0[:, jj], max_itr, tol)
            except (np.linalg.LinAlgError, ValueError, FloatingPointError):
                continue
            if not np.all(np.isfinite(gu)):
                continue
            g = _sign_unit(gu)
            try:
                ob = _cap_dir_obj(g, Sigma, H, X, Tvec, cov_shrinkage, B0, max_itr, tol)
            except RuntimeError:
                continue
            cands.append(g)
            objs.append(ob)
    if not cands:
        raise RuntimeError("cap_D1_opt: no initial value of gamma gave an estimate.")
    # the candidates in increasing order of the objective, the first one whose final
    # beta can be estimated (the best one when all succeed)
    for idx in np.argsort(np.asarray(objs), kind="stable"):
        gamma = cands[idx]
        sc, gg = scores(Sigma, gamma), float(gamma @ gamma)
        try:
            fin = _cap_beta_fit(sc, gg, X, Tvec, cov_shrinkage, B0, max_itr, tol)
        except RuntimeError:
            continue
        out = dict(gamma=gamma, **_beta_result(sc, gg, X, Tvec, cov_shrinkage, fin))
        return out, H
    raise RuntimeError("cap_D1_opt: beta could not be estimated for any initial value of gamma.")


def _cap_dk(Y, X, Phi0, cov_shrinkage, ninitial=None, seed=500, max_itr=1000, tol=1e-4,
            beta0=None, gamma0_mat=None):
    """Estimate the next direction after projecting out columns of ``Phi0``."""
    Y = as_matrix_list(Y)
    p = Y[0].shape[1]
    p0 = Phi0.shape[1]
    Sigma, Tvec = sigma_cube(Y, cov_shrinkage=cov_shrinkage)
    B0 = _beta0_mat(X.shape[1], beta0)
    # intercepts of the previous directions (beta fits from every initial value)
    beta_int0 = np.array([
        _cap_beta_fit(scores(Sigma, Phi0[:, j]), float(Phi0[:, j] @ Phi0[:, j]), X, Tvec,
                      cov_shrinkage, B0, max_itr, tol)["beta"][0]
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
            dnew = np.concatenate([d[: (p - p0)], np.sqrt(np.exp(beta_int0) * ntmp)])
            Ytmp.append((U[:, : len(dnew)] * dnew) @ Vt[: len(dnew), :])
    res, _ = _cap_d1_opt(Ytmp, X, cov_shrinkage, ninitial, seed, max_itr, tol,
                         beta0=beta0, gamma0_mat=gamma0_mat)
    res["orthogonal"] = Phi0.T @ res["gamma"]
    return res


# ------------------------------------------------------------- main entry -----
def cap_reg(Y, X, stop_crt="nD", nD=None, DfD_thred=5.0, cov_shrinkage=False,
            ninitial=None, seed=500, max_itr=1000, tol=1e-4, verbose=False,
            beta0=None, gamma0_mat=None):
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
    ninitial : number of initial values of gamma (default ``min(p + 6, 10)``).
    beta0 : initial value(s) of beta: a vector of length q, or a ``q x m`` array
            (one per column); ``None`` = two, ``0`` and :func:`cap_beta0_v5`.
            ``rep(0, q)`` alone was the start of capcov 0.1.0.
    gamma0_mat : ``p x k`` initial values of gamma (default: random, numpy).

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

    re1, _ = _cap_d1_opt(Y, X, cov_shrinkage, ninitial, seed, max_itr, tol,
                         beta0=beta0, gamma0_mat=gamma0_mat)
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
                rk = _cap_dk(Y, X, Phi, cov_shrinkage, ninitial, seed, max_itr, tol,
                             beta0=beta0, gamma0_mat=gamma0_mat)
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
                  seed=100, max_itr=1000, tol=1e-4, beta0=None):
    """Subject-level bootstrap inference for ``beta`` at a fixed ``gamma``.

    ``beta0``: initial value(s) of beta of every bootstrap fit (default: the two
    of :func:`cap_reg`); each fit keeps the converged estimate with the lowest
    objective.
    """
    Y = as_matrix_list(Y)
    X = np.asarray(X, float)
    n, q = len(Y), X.shape[1]
    B0 = _beta0_mat(q, beta0)
    boot = np.full((sims, q), np.nan)
    for b in range(sims):
        rng = np.random.RandomState(seed + b)
        idx = rng.randint(0, n, size=n)
        Yb = [Y[i] for i in idx]
        Xb = X[idx, :]
        try:
            res = cap_beta(Yb, Xb, gamma, cov_shrinkage=cov_shrinkage,
                           max_itr=max_itr, tol=tol, beta0=B0)
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
