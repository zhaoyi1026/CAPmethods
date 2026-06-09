"""CAP-CoC -- Covariance-on-Covariance regression (Python port of CAP-CoC/V3).

Model: for a pair of covariance-valued objects (outcome Y_i, predictor X_i) and
covariates W_i,

    log(gamma' Sy_i gamma) = alpha * log(theta' Sx_i theta) + W_i' beta ,

where ``gamma`` (length q) projects the outcome covariance and ``theta``
(length p) projects the predictor covariance. ``coc_reg`` finds the leading k
projection pairs (gamma_k, theta_k) with per-direction (alpha_k, beta_k).

Ported from ``CAP-CoC/V3/COCReg.R`` (RcppArmadillo-accelerated R reference);
pure numpy/scipy here. The R version remains the reference.

Public API
----------
``coc_reg(Y, X, W, ...)``           -- find directions (== R ``COCReg``)
``coc_coef(Y, X, W, gamma, theta)`` -- alpha/beta given gamma, theta (== ``COCReg.coef``)
``coc_coef_asmp(...)``              -- asymptotic inference (== ``COCReg.coef.asmp``)
``coc_coef_boot(...)``              -- bootstrap inference (== ``COCReg.coef.boot``)
"""
from __future__ import annotations

import numpy as np

from ._core import as_matrix_list, cov_ls, scores, accum, gamma_solve, diag_level

__all__ = ["coc_reg", "coc_coef", "coc_coef_asmp", "coc_coef_boot"]


# --------------------------------------------------------------------------- #
# covariance cubes (shape (p, p, n)), matching coc_kernels.cpp                 #
# --------------------------------------------------------------------------- #
def _sample_cov_cube(data):
    """Per-subject centered sample covariance ``X_i' X_i / T_i`` -> ``(p, p, n)``."""
    data = as_matrix_list(data)
    n = len(data)
    p = data[0].shape[1]
    cube = np.empty((p, p, n))
    for i, Xi in enumerate(data):
        Xc = Xi - Xi.mean(axis=0, keepdims=True)
        cube[:, :, i] = (Xc.T @ Xc) / Xi.shape[0]
    return cube


def _covls_cube(Y):
    """Per-subject linear-shrinkage (cov.ls) covariance -> ``(q, q, n)``."""
    Y = as_matrix_list(Y)
    n = len(Y)
    q = Y[0].shape[1]
    cube = np.empty((q, q, n))
    for i, Yi in enumerate(Y):
        cube[:, :, i] = cov_ls(Yi)
    return cube


def _nf2(A):
    """Standardized squared Frobenius norm ``sum(A_ij^2) / nrow(A)``."""
    return np.sum(A * A) / A.shape[0]


def _cov_sk_x(X):
    """Shrinkage covariance cube for the predictor X (== R ``cov.sk.x``)."""
    X = as_matrix_list(X)
    n = len(X)
    p = X[0].shape[1]
    S = np.empty((p, p, n))
    nx = np.empty(n)
    nuv = np.empty(n)
    Xc = []
    for i, Xi in enumerate(X):
        nx[i] = Xi.shape[0]
        xc = Xi - Xi.mean(axis=0, keepdims=True)
        Xc.append(xc)
        Si = (xc.T @ xc) / nx[i]
        S[:, :, i] = Si
        nuv[i] = np.trace(Si) / p
    nuhat = nuv.mean()
    tau2 = np.empty(n)
    omega2 = np.empty(n)
    eps2 = np.empty(n)
    Ip = np.eye(p)
    for i in range(n):
        Si = S[:, :, i]
        tau2[i] = _nf2(Si - nuhat * Ip)
        normS2 = np.sum(Si * Si)
        xx = np.einsum("tj,tj->t", Xc[i], Xc[i])
        xSx = np.einsum("tj,jk,tk->t", Xc[i], Si, Xc[i])
        otmp = np.sum((xx * xx - 2.0 * xSx + normS2) / p) / (nx[i] * nx[i])
        omega2[i] = min(otmp, tau2[i])
        eps2[i] = tau2[i] - omega2[i]
    tauh, omh, eph = tau2.mean(), omega2.mean(), eps2.mean()
    out = np.empty((p, p, n))
    for i in range(n):
        out[:, :, i] = (omh / tauh) * nuhat * Ip + (eph / tauh) * S[:, :, i]
    return out


def _cov_sk_y(Y, gamma, kappa):
    """Shrinkage covariance cube for the outcome Y given gamma, kappa (== ``cov.sk.y``)."""
    Y = as_matrix_list(Y)
    gamma = np.asarray(gamma, dtype=float).ravel()
    kappa = np.asarray(kappa, dtype=float).ravel()
    n = len(Y)
    q = Y[0].shape[1]
    Sy = np.empty((q, q, n))
    ny = np.empty(n)
    for i, Yi in enumerate(Y):
        ny[i] = Yi.shape[0]
        yc = Yi - Yi.mean(axis=0, keepdims=True)
        Sy[:, :, i] = (yc.T @ yc) / ny[i]
    gg = float(gamma @ gamma)
    muhat = kappa.mean() / gg
    d2 = np.empty(n)
    psi2 = np.empty(n)
    phi2 = np.empty(n)
    for i in range(n):
        gs = float(gamma @ Sy[:, :, i] @ gamma)
        d2[i] = (gs - muhat * gg) ** 2
        psi2[i] = min((gs - kappa[i]) ** 2 / ny[i], d2[i])
        phi2[i] = d2[i] - psi2[i]
    dh, ph, ps = d2.mean(), phi2.mean(), psi2.mean()
    out = np.empty((q, q, n))
    Iq = np.eye(q)
    for i in range(n):
        out[:, :, i] = (ps / dh) * muhat * Iq + (ph / dh) * Sy[:, :, i]
    return out


# --------------------------------------------------------------------------- #
# small helpers                                                                #
# --------------------------------------------------------------------------- #
def _as_W(W):
    W = np.asarray(W, dtype=float)
    if W.ndim == 1:
        W = W[:, None]
    return W


def _sizes(data):
    return np.array([np.asarray(d).shape[0] for d in data], dtype=float)


def _Sx_cube(X, cov_shrinkage_x):
    return _cov_sk_x(X) if cov_shrinkage_x else _sample_cov_cube(X)


def _Sy_cube(Y, cov_shrinkage_y):
    return _covls_cube(Y) if cov_shrinkage_y else _sample_cov_cube(Y)


def _obj_func(Sy, Sx, W, gamma, theta, alpha, beta):
    sy = scores(Sy, gamma)
    sx = scores(Sx, theta)
    lv = (np.log(sy) - alpha * np.log(sx) - W @ beta) ** 2
    return np.nanmean(lv)


def _sign_norm(v):
    v = v / np.sqrt(np.sum(v ** 2))
    if v[np.argmax(np.abs(v))] < 0:
        v = -v
    return v


# --------------------------------------------------------------------------- #
# given gamma, theta: estimate alpha, beta                                     #
# --------------------------------------------------------------------------- #
def coc_coef(Y, X, W, gamma, theta, cov_shrinkage_y=True, cov_shrinkage_x=True,
             max_itr=1000, tol=1e-4):
    """Estimate (alpha, beta) given fixed (gamma, theta) (== R ``COCReg.coef``)."""
    Y = as_matrix_list(Y)
    X = as_matrix_list(X)
    W = _as_W(W)
    n = len(Y)
    p = X[0].shape[1]
    q = Y[0].shape[1]
    r = W.shape[1]
    gamma = np.asarray(gamma, dtype=float).ravel()
    theta = np.asarray(theta, dtype=float).ravel()

    nxvec = _sizes(X)
    nyvec = _sizes(Y)
    if nxvec.min() - 5 < p:
        cov_shrinkage_x = True
    if nyvec.min() - 5 < q:
        cov_shrinkage_y = True

    Sx = _Sx_cube(X, cov_shrinkage_x)
    Sy = _Sy_cube(Y, cov_shrinkage_y)
    score_x = scores(Sx, theta)
    score_y = scores(Sy, gamma)

    alpha0 = 0.0
    beta0 = np.zeros(r)          # init irrelevant to the converged estimate
    WtW_inv = np.linalg.pinv(W.T @ W / n)

    s = 0
    diff = 100.0
    alpha_new, beta_new = alpha0, beta0
    while s <= max_itr and diff > tol:
        s += 1
        lsx = np.log(score_x)
        lsy = np.log(score_y)
        alpha_new = np.mean(lsy * lsx - (W @ beta0) * lsx) / np.mean(lsx ** 2)
        otmp = ((lsy - alpha_new * lsx)[:, None] * W).mean(axis=0)
        beta_new = WtW_inv @ otmp
        if cov_shrinkage_y:
            kappa = np.exp(alpha_new * lsx + W @ beta_new)
            Sy = _cov_sk_y(Y, gamma, kappa)
            score_y = scores(Sy, gamma)
        diff = max(np.max(np.abs(beta_new - beta0)), abs(alpha_new - alpha0))
        beta0, alpha0 = beta_new, alpha_new

    return {
        "gamma": gamma, "theta": theta, "alpha": float(alpha_new),
        "beta": beta_new, "convergence": s < max_itr,
        "score_y": score_y, "score_x": score_x,
    }


# --------------------------------------------------------------------------- #
# one direction: estimate gamma, theta, alpha, beta                           #
# --------------------------------------------------------------------------- #
def _coc_d1_base(Y, X, W, gamma=None, theta=None, Hy=None, Hx=None,
                 cov_shrinkage_y=True, cov_shrinkage_x=True,
                 max_itr=1000, tol=1e-4, gamma0=None, theta0=None, burn_in=1000):
    Y = as_matrix_list(Y)
    X = as_matrix_list(X)
    W = _as_W(W)
    n = len(Y)
    p = X[0].shape[1]
    q = Y[0].shape[1]
    r = W.shape[1]

    nxvec = _sizes(X)
    nyvec = _sizes(Y)
    if nxvec.min() - 5 < p:
        cov_shrinkage_x = True
    if nyvec.min() - 5 < q:
        cov_shrinkage_y = True

    Sx = _Sx_cube(X, cov_shrinkage_x)
    Sy = _Sy_cube(Y, cov_shrinkage_y)

    if Hx is None:
        Hx = accum(Sx, nxvec / nxvec.sum())
    if Hy is None:
        Hy = accum(Sy, nyvec / nyvec.sum())

    alpha0 = 0.0
    beta0 = np.zeros(r)
    if gamma0 is None and gamma is None:
        gamma0 = np.full(q, 1.0 / np.sqrt(q))
    elif gamma0 is None:
        gamma0 = np.asarray(gamma, dtype=float).ravel()
    else:
        gamma0 = np.asarray(gamma0, dtype=float).ravel()
    if theta0 is None and theta is None:
        theta0 = np.full(p, 1.0 / np.sqrt(p))
    elif theta0 is None:
        theta0 = np.asarray(theta, dtype=float).ravel()
    else:
        theta0 = np.asarray(theta0, dtype=float).ravel()

    gamma_fixed = None if gamma is None else np.asarray(gamma, dtype=float).ravel()
    theta_fixed = None if theta is None else np.asarray(theta, dtype=float).ravel()

    score_x = scores(Sx, theta0)
    score_y = scores(Sy, gamma0)
    WtW_inv = np.linalg.pinv(W.T @ W / n)

    obj0 = _obj_func(Sy, Sx, W, gamma0, theta0, alpha0, beta0)

    s = 0
    diff = 100.0
    while s <= max_itr and diff > tol:
        s += 1
        lsx = np.log(score_x)
        lsy = np.log(score_y)

        alpha_new = np.mean(lsy * lsx - (W @ beta0) * lsx) / np.mean(lsx ** 2)
        otmp = ((lsy - alpha_new * lsx)[:, None] * W).mean(axis=0)
        beta_new = WtW_inv @ otmp

        if theta_fixed is None:
            Vtmp = lsy - W @ beta_new
            w2 = ((alpha_new * lsx - Vtmp) / score_x) * (2.0 * alpha_new / n)
            A2 = accum(Sx, w2)
            theta_new = gamma_solve(A2, Hx)
            score_x = scores(Sx, theta_new)
        else:
            theta_new = theta_fixed

        if gamma_fixed is None:
            lsx_u = np.log(score_x)               # theta may have changed
            Utmp = alpha_new * lsx_u + W @ beta_new
            if cov_shrinkage_y:
                kappa = np.exp(alpha_new * lsx_u + W @ beta_new)
                Sy = _cov_sk_y(Y, gamma0, kappa)
            w1 = ((np.log(score_y) - Utmp) / score_y) * (2.0 / n)
            A1 = accum(Sy, w1)
            gamma_new = gamma_solve(A1, Hy)
            score_y = scores(Sy, gamma_new)
        else:
            gamma_new = gamma_fixed

        diff1 = max(np.max(np.abs(gamma_new - gamma0)), np.max(np.abs(theta_new - theta0)))
        diff2 = max(np.max(np.abs(beta_new - beta0)), abs(alpha_new - alpha0))
        diff = max(diff1, diff2)

        obj_new = _obj_func(Sy, Sx, W, gamma_new, theta_new, alpha_new, beta_new)
        if obj_new > obj0 and s > burn_in:
            break
        obj0 = obj_new
        beta0, alpha0, gamma0, theta0 = beta_new, alpha_new, gamma_new, theta_new

    gamma_est = _sign_norm(gamma0)
    theta_est = _sign_norm(theta0)
    otmp = coc_coef(Y, X, W, gamma_est, theta_est, cov_shrinkage_y=cov_shrinkage_y,
                    cov_shrinkage_x=cov_shrinkage_x, max_itr=max_itr, tol=tol)
    return {
        "gamma": gamma_est, "theta": theta_est, "alpha": otmp["alpha"],
        "beta": otmp["beta"], "convergence": s < max_itr,
        "score_y": otmp["score_y"], "score_x": otmp["score_x"],
    }


def _coc_d1(Y, X, W, gamma=None, theta=None, Hy=None, Hx=None,
            cov_shrinkage_y=True, cov_shrinkage_x=True, max_itr=1000, tol=1e-4,
            gamma0=None, theta0=None, burn_in=1000):
    if gamma is None and theta is None:
        tmp = _coc_d1_base(Y, X, W, gamma=None, theta=None, Hy=Hy, Hx=Hx,
                           cov_shrinkage_y=cov_shrinkage_y, cov_shrinkage_x=cov_shrinkage_x,
                           max_itr=max_itr, tol=tol, gamma0=gamma0, theta0=theta0, burn_in=burn_in)
        return _coc_d1_base(Y, X, W, gamma=tmp["gamma"], theta=None, Hy=Hy, Hx=Hx,
                            cov_shrinkage_y=cov_shrinkage_y, cov_shrinkage_x=cov_shrinkage_x,
                            max_itr=max_itr, tol=tol, gamma0=gamma0, theta0=theta0, burn_in=burn_in)
    return _coc_d1_base(Y, X, W, gamma=gamma, theta=theta, Hy=Hy, Hx=Hx,
                        cov_shrinkage_y=cov_shrinkage_y, cov_shrinkage_x=cov_shrinkage_x,
                        max_itr=max_itr, tol=tol, gamma0=gamma0, theta0=theta0, burn_in=burn_in)


def _eig_desc(S):
    """Eigenvectors of symmetric ``S`` in descending eigenvalue order (R ``eigen``)."""
    w, V = np.linalg.eigh(S)
    return V[:, ::-1]


def _coc_d1_opt(Y, X, W, gamma=None, theta=None, Hy=None, Hx=None,
                cov_shrinkage_y=True, cov_shrinkage_x=True, max_itr=1000, tol=1e-4,
                burn_in=1000, gamma0_mat=None, theta0_mat=None, ninitial=None):
    Y = as_matrix_list(Y)
    X = as_matrix_list(X)
    W = _as_W(W)
    n = len(Y)
    p = X[0].shape[1]
    q = Y[0].shape[1]

    nxvec = _sizes(X)
    nyvec = _sizes(Y)
    if nxvec.min() - 5 < p:
        cov_shrinkage_x = True
    if nyvec.min() - 5 < q:
        cov_shrinkage_y = True

    Sx = _Sx_cube(X, cov_shrinkage_x)
    Sy = _Sy_cube(Y, cov_shrinkage_y)
    Sx_bar = Sx.mean(axis=2)
    Sy_bar = Sy.mean(axis=2)

    if Hx is None:
        Hx = accum(Sx, nxvec / nxvec.sum())
    if Hy is None:
        Hy = accum(Sy, nyvec / nyvec.sum())

    if ninitial is None:
        if gamma0_mat is not None and theta0_mat is not None:
            ninitial = min(gamma0_mat.shape[1], theta0_mat.shape[1], 10)
        elif gamma0_mat is not None:
            ninitial = min(gamma0_mat.shape[1], 10)
        elif theta0_mat is not None:
            ninitial = min(theta0_mat.shape[1], 10)
        else:
            ninitial = min(p, q, 10)

    if gamma is None:
        if gamma0_mat is None:
            gamma0_mat = _eig_desc(Sy_bar)[:, :ninitial]
    else:
        g = np.asarray(gamma, dtype=float).ravel()
        gamma0_mat = np.tile(g[:, None], (1, ninitial))
    if theta is None:
        if theta0_mat is None:
            theta0_mat = _eig_desc(Sx_bar)[:, :ninitial]
    else:
        t = np.asarray(theta, dtype=float).ravel()
        theta0_mat = np.tile(t[:, None], (1, ninitial))

    best = None
    best_obj = np.inf
    for kk in range(ninitial):
        try:
            rk = _coc_d1(Y, X, W, gamma=gamma, theta=theta, Hy=Hy, Hx=Hx,
                         cov_shrinkage_y=cov_shrinkage_y, cov_shrinkage_x=cov_shrinkage_x,
                         max_itr=max_itr, tol=tol, gamma0=gamma0_mat[:, kk],
                         theta0=theta0_mat[:, kk], burn_in=burn_in)
        except Exception:
            continue
        g = rk["gamma"]
        t = rk["theta"]
        g_un = g / np.sqrt(g @ Hy @ g)
        t_un = t / np.sqrt(t @ Hx @ t)
        try:
            coef = coc_coef(Y, X, W, g_un, t_un, cov_shrinkage_y=cov_shrinkage_y,
                            cov_shrinkage_x=cov_shrinkage_x, max_itr=max_itr, tol=tol)
            obj = _obj_func(Sy, Sx, W, g_un, t_un, coef["alpha"], coef["beta"])
        except Exception:
            obj = np.inf
        if obj < best_obj:
            best_obj = obj
            best = rk
    return best


def _coc_dk(Y, X, W, Gamma0=None, Theta0=None, remove_y=True, remove_x=True,
            Hy=None, Hx=None, cov_shrinkage_y=True, cov_shrinkage_x=True,
            max_itr=1000, tol=1e-4, burn_in=1000, ninitial=None):
    Y = as_matrix_list(Y)
    X = as_matrix_list(X)
    if Gamma0 is None and Theta0 is None:
        return _coc_d1_opt(Y, X, W, Hy=Hy, Hx=Hx, cov_shrinkage_y=cov_shrinkage_y,
                           cov_shrinkage_x=cov_shrinkage_x, max_itr=max_itr, tol=tol,
                           burn_in=burn_in, ninitial=ninitial)
    if Gamma0 is not None and remove_y:
        Py = Gamma0 @ Gamma0.T
        Ytmp = [Yi - Yi @ Py for Yi in Y]
    else:
        Ytmp = Y
    if Theta0 is not None and remove_x:
        Px = Theta0 @ Theta0.T
        Xtmp = [Xi - Xi @ Px for Xi in X]
    else:
        Xtmp = X
    re = _coc_d1_opt(Ytmp, Xtmp, W, Hy=Hy, Hx=Hx, cov_shrinkage_y=cov_shrinkage_y,
                     cov_shrinkage_x=cov_shrinkage_x, max_itr=max_itr, tol=tol,
                     burn_in=burn_in, ninitial=ninitial)
    if re is not None and Gamma0 is not None:
        re["orthogonal_gamma"] = re["gamma"] @ Gamma0
        re["orthogonal_theta"] = re["theta"] @ Theta0
    return re


# --------------------------------------------------------------------------- #
# main entry: find the first k directions                                      #
# --------------------------------------------------------------------------- #
def coc_reg(Y, X, W, stop_crt="nD", nD=None, DfD_y_thred=2.0, DfD_thred=2.0,
            remove_y=True, remove_x=True, Hy=None, Hx=None,
            cov_shrinkage_y=True, cov_shrinkage_x=True, max_itr=1000, tol=1e-4,
            burn_in=1000, ninitial=None, verbose=False):
    """Covariance-on-covariance regression (== R ``COCReg``).

    Parameters
    ----------
    Y, X : list of arrays
        ``Y[i]`` is ``Ty_i x q``; ``X[i]`` is ``Tx_i x p``.
    W : (n, r) array
        Covariates (include an intercept column).
    stop_crt : {"nD", "DfD.y", "DfD"}
        Number of directions, or a deviation-from-diagonality threshold.

    Returns a dict with ``gamma`` (q x k), ``theta`` (p x k), ``alpha`` (k,),
    ``beta`` (r x k), ``score_y``/``score_x`` (n x k), and DfD diagnostics.
    """
    Y = as_matrix_list(Y)
    X = as_matrix_list(X)
    W = _as_W(W)
    n = len(Y)
    p = X[0].shape[1]
    q = Y[0].shape[1]

    if stop_crt == "nD" and nD is None:
        stop_crt = "DfD"

    nxvec = _sizes(X)
    nyvec = _sizes(Y)
    if nxvec.min() - 5 < p:
        cov_shrinkage_x = True
    if nyvec.min() - 5 < q:
        cov_shrinkage_y = True

    re1 = _coc_d1_opt(Y, X, W, Hy=Hy, Hx=Hx, cov_shrinkage_y=cov_shrinkage_y,
                      cov_shrinkage_x=cov_shrinkage_x, max_itr=max_itr, tol=tol,
                      burn_in=burn_in, ninitial=ninitial)
    Gamma = [re1["gamma"]]
    Theta = [re1["theta"]]
    alpha = [re1["alpha"]]
    beta = [re1["beta"]]
    score_y = [re1["score_y"]]
    score_x = [re1["score_x"]]
    if verbose:
        print("Component 1")

    def _add(re_tmp):
        Gamma.append(re_tmp["gamma"])
        Theta.append(re_tmp["theta"])
        alpha.append(re_tmp["alpha"])
        beta.append(re_tmp["beta"])
        score_y.append(re_tmp["score_y"])
        score_x.append(re_tmp["score_x"])
        if verbose:
            print(f"Component {len(Gamma)}")

    if stop_crt == "nD":
        for _ in range(2, nD + 1):
            re_tmp = _coc_dk(Y, X, W, Gamma0=np.column_stack(Gamma), Theta0=np.column_stack(Theta),
                             remove_y=remove_y, remove_x=remove_x, Hy=Hy, Hx=Hx,
                             cov_shrinkage_y=cov_shrinkage_y, cov_shrinkage_x=cov_shrinkage_x,
                             max_itr=max_itr, tol=tol, burn_in=burn_in, ninitial=ninitial)
            if re_tmp is None:
                break
            _add(re_tmp)
    elif stop_crt == "DfD.y":
        DfD_tmp = 1.0
        while DfD_tmp < DfD_y_thred:
            re_tmp = _coc_dk(Y, X, W, Gamma0=np.column_stack(Gamma), Theta0=np.column_stack(Theta),
                             remove_y=remove_y, remove_x=remove_x, Hy=Hy, Hx=Hx,
                             cov_shrinkage_y=cov_shrinkage_y, cov_shrinkage_x=cov_shrinkage_x,
                             max_itr=max_itr, tol=tol, burn_in=burn_in, ninitial=ninitial)
            if re_tmp is None:
                break
            k = len(Gamma)   # 0-based index of the candidate new direction
            avg_y, _ = diag_level(Y, np.column_stack(Gamma + [re_tmp["gamma"]]))
            DfD_tmp = avg_y[k]
            if DfD_tmp < DfD_y_thred:
                _add(re_tmp)
            else:
                break
    elif stop_crt == "DfD":
        DfD_tmp = 1.0
        while DfD_tmp < DfD_thred:
            re_tmp = _coc_dk(Y, X, W, Gamma0=np.column_stack(Gamma), Theta0=np.column_stack(Theta),
                             remove_y=remove_y, remove_x=remove_x, Hy=Hy, Hx=Hx,
                             cov_shrinkage_y=cov_shrinkage_y, cov_shrinkage_x=cov_shrinkage_x,
                             max_itr=max_itr, tol=tol, burn_in=burn_in, ninitial=ninitial)
            if re_tmp is None:
                break
            k = len(Gamma)
            avg_y, _ = diag_level(Y, np.column_stack(Gamma + [re_tmp["gamma"]]))
            avg_x, _ = diag_level(X, np.column_stack(Theta + [re_tmp["theta"]]))
            DfD_tmp = np.nanmax([avg_y[k], avg_x[k]])
            if DfD_tmp < DfD_thred and not np.any(np.isnan(np.r_[avg_y, avg_x])):
                _add(re_tmp)
            else:
                break
    else:
        raise ValueError("stop_crt must be 'nD', 'DfD.y', or 'DfD'")

    Gmat = np.column_stack(Gamma)
    Tmat = np.column_stack(Theta)
    avg_y, sub_y = diag_level(Y, Gmat)
    avg_x, sub_x = diag_level(X, Tmat)
    cols = [f"C{i + 1}" for i in range(Gmat.shape[1])]

    return {
        "gamma": Gmat, "theta": Tmat,
        "alpha": np.array(alpha), "beta": np.column_stack(beta),
        "score_y": np.column_stack(score_y), "score_x": np.column_stack(score_x),
        "DfD_y": {"avg_level": avg_y, "sub_level": sub_y},
        "DfD_x": {"avg_level": avg_x, "sub_level": sub_x},
        "orthogonality_y": Gmat.T @ Gmat, "orthogonality_x": Tmat.T @ Tmat,
        "columns": cols,
    }


# --------------------------------------------------------------------------- #
# inference                                                                    #
# --------------------------------------------------------------------------- #
def coc_coef_asmp(Y, X, W, gamma, theta, conf_level=0.95, cov_shrinkage_y=True,
                  cov_shrinkage_x=True, max_itr=1000, tol=1e-4):
    """Asymptotic inference for (alpha, beta) (== R ``COCReg.coef.asmp``)."""
    from scipy.stats import norm
    Y = as_matrix_list(Y)
    X = as_matrix_list(X)
    W = _as_W(W)
    n = len(Y)
    p = X[0].shape[1]
    q = Y[0].shape[1]
    nyvec = _sizes(Y)
    nxvec = _sizes(X)
    if nxvec.min() - 5 < p:
        cov_shrinkage_x = True
    if nyvec.min() - 5 < q:
        cov_shrinkage_y = True

    Sx = _Sx_cube(X, cov_shrinkage_x)
    theta = np.asarray(theta, dtype=float).ravel()
    gamma = np.asarray(gamma, dtype=float).ravel()
    score_x = scores(Sx, theta)

    coef = coc_coef(Y, X, W, gamma, theta, cov_shrinkage_y=cov_shrinkage_y,
                    cov_shrinkage_x=cov_shrinkage_x, max_itr=max_itr, tol=tol)

    lsx = np.log(score_x)
    Gx = np.nanmean(lsx ** 2)
    Qw = W.T @ W / n
    Hxw = np.nanmean(lsx[:, None] * W, axis=0)
    Fs = np.block([[np.array([[Gx]]), Hxw[None, :]],
                   [Hxw[:, None], Qw]])
    cov_asmp = np.linalg.pinv(Fs) / nyvec.sum()
    se = np.sqrt(np.diag(cov_asmp))
    z = norm.ppf(1 - (1 - conf_level) / 2)

    a_est = coef["alpha"]
    a_se = se[0]
    a_stat = a_est / a_se
    a_pv = (1 - norm.cdf(abs(a_stat))) * 2
    b_est = coef["beta"]
    b_se = se[1:]
    b_stat = b_est / b_se
    b_pv = (1 - norm.cdf(np.abs(b_stat))) * 2

    return {
        "alpha": {"estimate": a_est, "se": a_se, "statistic": a_stat,
                  "pvalue": a_pv, "lower": a_est - z * a_se, "upper": a_est + z * a_se},
        "beta": {"estimate": b_est, "se": b_se, "statistic": b_stat,
                 "pvalue": b_pv, "lower": b_est - z * b_se, "upper": b_est + z * b_se},
        "cov_asmp": cov_asmp,
    }


def coc_coef_boot(Y, X, W, gamma, theta, cov_shrinkage_y=True, cov_shrinkage_x=True,
                  sims=1000, boot_ci_type="se", conf_level=0.95, max_itr=1000,
                  tol=1e-4, seed=100):
    """Bootstrap inference for (alpha, beta) (== R ``COCReg.coef.boot``).

    Note: R seeds each resample with ``set.seed(100+b)``; this uses a numpy
    ``RandomState`` so the specific resamples differ, but the inference is
    statistically equivalent.
    """
    from scipy.stats import norm
    Y = as_matrix_list(Y)
    X = as_matrix_list(X)
    W = _as_W(W)
    n = len(Y)
    r = W.shape[1]
    gamma = np.asarray(gamma, dtype=float).ravel()
    theta = np.asarray(theta, dtype=float).ravel()

    alpha_boot = np.full(sims, np.nan)
    beta_boot = np.full((r, sims), np.nan)
    for b in range(sims):
        rng = np.random.RandomState(seed + b)
        idx = rng.randint(0, n, size=n)
        Yt = [Y[i] for i in idx]
        Xt = [X[i] for i in idx]
        Wt = W[idx, :]
        try:
            rt = coc_coef(Yt, Xt, Wt, gamma, theta, cov_shrinkage_y=cov_shrinkage_y,
                          cov_shrinkage_x=cov_shrinkage_x, max_itr=max_itr, tol=tol)
            alpha_boot[b] = rt["alpha"]
            beta_boot[:, b] = rt["beta"]
        except Exception:
            pass

    a_est = np.nanmean(alpha_boot)
    a_se = np.nanstd(alpha_boot, ddof=1)
    b_est = np.nanmean(beta_boot, axis=1)
    b_se = np.nanstd(beta_boot, axis=1, ddof=1)
    a_stat = a_est / a_se
    b_stat = b_est / b_se
    a_pv = (1 - norm.cdf(abs(a_stat))) * 2
    b_pv = (1 - norm.cdf(np.abs(b_stat))) * 2

    if boot_ci_type == "se":
        z = norm.ppf(1 - (1 - conf_level) / 2)
        a_ci = (a_est - z * a_se, a_est + z * a_se)
        b_ci = np.column_stack([b_est - z * b_se, b_est + z * b_se])
    else:
        lo, hi = (1 - conf_level) / 2, 1 - (1 - conf_level) / 2
        a_ci = tuple(np.nanquantile(alpha_boot, [lo, hi]))
        b_ci = np.column_stack([np.nanquantile(beta_boot, lo, axis=1),
                                np.nanquantile(beta_boot, hi, axis=1)])

    return {
        "alpha": {"estimate": a_est, "se": a_se, "statistic": a_stat,
                  "pvalue": a_pv, "lower": a_ci[0], "upper": a_ci[1]},
        "beta": {"estimate": b_est, "se": b_se, "statistic": b_stat,
                 "pvalue": b_pv, "lower": b_ci[:, 0], "upper": b_ci[:, 1]},
        "alpha_boot": alpha_boot, "beta_boot": beta_boot,
    }
