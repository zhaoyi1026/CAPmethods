"""LCAP (gamma-invariant longitudinal CAP) -- Python port of LCAP_gamma-invar/V6.

Longitudinal CAP where a single projection ``gamma`` is shared across visits and
subjects, and the projected log-variance follows a subject-level random-effects
(random intercept + random slopes) model:

    log(gamma' Sigma_{iv} gamma) = x_{iv}' b_i ,
    b_i = (beta0_i, beta1_i) ~ N(beta, diag(sigma2, Omega)) .

``cap_reg`` finds the leading directions; ``cap_beta`` estimates the random-
effects model for a fixed gamma; ``cap_beta_boot`` does bootstrap inference.

Ported from ``LCAP_gamma-invar/V6/Longitudinal_HDCAP.R`` and its RcppArmadillo
kernel ``cap_kernels.cpp`` (the recursion math is taken verbatim from the
kernel). The R version remains the reference. ``cov_beta_diag=True`` (diagonal
random-effect covariance) is the supported mode, matching the web app.

Data layout
-----------
``Y[i][v]`` is a ``T_{iv} x p`` matrix (subject ``i``, visit ``v``);
``X[i]`` is an ``nV_i x q`` design matrix whose first column is the intercept.
"""
from __future__ import annotations

import numpy as np

from ._core import cov_ls, gamma_solve

__all__ = ["cap_beta", "cap_reg", "cap_beta_boot", "diag_level"]


# --------------------------------------------------------------------------- #
# low-level helpers (mirror cap_kernels.cpp)                                   #
# --------------------------------------------------------------------------- #
def _ginv_diag(d):
    """MASS::ginv of a diagonal PSD matrix: reciprocal-or-zero."""
    d = np.asarray(d, dtype=float)
    r = np.zeros_like(d)
    if d.size:
        tol = np.sqrt(np.finfo(float).eps) * d.max()
        m = d > tol
        r[m] = 1.0 / d[m]
    return r


def _visit_covs(Y, shrinkage):
    """Per-subject per-visit covariances and visit sample sizes.

    Returns ``S`` (list of ``(p, p, nV_i)`` arrays) and ``nT`` (list of length-nV_i
    arrays). Sample covariance is ``cov(Y_iv) * (T-1)/T`` (centered, divided by T).
    """
    n = len(Y)
    p = Y[0][0].shape[1]
    S = []
    nT = []
    for i in range(n):
        nv = len(Y[i])
        Si = np.empty((p, p, nv))
        ni = np.empty(nv)
        for v in range(nv):
            Yiv = np.asarray(Y[i][v], dtype=float)
            T = Yiv.shape[0]
            ni[v] = T
            if shrinkage:
                Si[:, :, v] = cov_ls(Yiv)
            else:
                Yc = Yiv - Yiv.mean(axis=0, keepdims=True)
                Si[:, :, v] = (Yc.T @ Yc) / T
        S.append(Si)
        nT.append(ni)
    return S, nT


def _scores_list(S, gamma):
    """Per-subject projected scores ``gamma' S_iv gamma`` as a list of vectors."""
    gamma = np.asarray(gamma, dtype=float).ravel()
    return [np.einsum("i,ijv,j->v", gamma, Si, gamma) for Si in S]


def _flatten(X, nT):
    """Stack the design and visit sizes across subjects -> (X_flat, Tvec, vcount)."""
    Xf = np.vstack([np.asarray(Xi, dtype=float) for Xi in X])
    Tflat = np.concatenate([np.asarray(t, dtype=float) for t in nT])
    vcount = np.array([np.asarray(Xi).shape[0] for Xi in X], dtype=int)
    return Xf, Tflat, vcount


def _init_re(scores, X, q):
    """Per-subject lm of ``log(score)`` on ``X[:, 1:]`` (== .lcap_initRE)."""
    n = len(X)
    beta0vec = np.empty(n)
    beta1mat = np.empty((n, q - 1))
    for i in range(n):
        Xi = np.asarray(X[i], dtype=float)
        y = np.log(scores[i])
        coef, *_ = np.linalg.lstsq(Xi, y, rcond=None)   # intercept is X[:,0]
        beta0vec[i] = coef[0]
        beta1mat[i, :] = coef[1:]
    beta0 = beta0vec.mean()
    beta1 = beta1mat.mean(axis=0)
    sigma2 = np.mean((beta0vec - beta0) ** 2)
    OmegaDiag = np.array([np.mean((beta1mat[:, j] - beta1[j]) ** 2)
                          for j in range(q - 1)])
    beta = np.concatenate([[beta0], beta1])
    return beta0vec, beta1mat, beta, sigma2, OmegaDiag


def _covlsconst_rho(scoreRaw, gg, Xf, Tvec, vcount, beta0vec, beta1mat, n, q):
    """rho1, rho2 of the constant-shrinkage covariance (== covlsconst_rho)."""
    # mu (with V5 nested averaging)
    mu = 0.0
    idx = 0
    for i in range(n):
        acc = 0.0
        bi = np.concatenate([[beta0vec[i]], beta1mat[i]])
        for _ in range(vcount[i]):
            acc += np.exp(Xf[idx] @ bi)
            idx += 1
        mu += acc
        mu /= vcount[i]
    mu /= (n * gg)
    hd = hp = hf = 0.0
    idx = 0
    for i in range(n):
        sd = sp = sf = 0.0
        bi = np.concatenate([[beta0vec[i]], beta1mat[i]])
        for _ in range(vcount[i]):
            u1 = Xf[idx] @ bi
            sc = scoreRaw[idx]
            dlt = (sc - mu * gg) ** 2
            psi = min((sc - np.exp(u1)) ** 2 / Tvec[idx], dlt)
            sd += dlt; sp += psi; sf += (dlt - psi)
            idx += 1
        hd += sd; hd /= vcount[i]
        hp += sp; hp /= vcount[i]
        hf += sf; hf /= vcount[i]
    hd /= n; hp /= n; hf /= n
    rho1 = hp * mu / hd
    rho2 = hf / hd
    return rho1, rho2, mu, hd, hp, hf


def _re_recursion(scoreWork, scoreRaw, gg, Xf, Tvec, vcount, shrink, n, q,
                  beta0vec, beta1mat, beta, sigma2, OmegaDiag, max_itr, tol):
    """Random-effects EM given fixed gamma (== re_recursion in the kernel)."""
    scoreWork = scoreWork.copy()
    beta0vec = beta0vec.copy(); beta1mat = beta1mat.copy(); beta = beta.copy()
    start = np.concatenate([[0], np.cumsum(vcount)[:-1]])
    s = 0
    diff = 100.0
    while s <= max_itr and diff > tol:
        s += 1
        dcov = np.concatenate([[sigma2], OmegaDiag])
        gdiag = _ginv_diag(dcov)
        bmat = np.empty((n, q))
        for i in range(n):
            b_i = np.concatenate([[beta0vec[i]], beta1mat[i]])
            pt1 = np.zeros(q)
            pt2 = np.zeros((q, q))
            for v in range(vcount[i]):
                r = start[i] + v
                xr = Xf[r]
                u1 = xr @ b_i
                e = np.exp(-u1)
                w = Tvec[r] / 2.0
                pt1 += (xr - scoreWork[r] * e * xr) * w
                pt2 += (scoreWork[r] * e) * np.outer(xr, xr) * w
            pt1 += gdiag * (b_i - beta)
            pt2[np.diag_indices(q)] += gdiag
            bmat[i] = b_i - np.linalg.solve(pt2, pt1)
        bmean = bmat.mean(axis=0)
        Bc = bmat - bmean
        bcov = (Bc.T @ Bc) / n
        beta_upd = bmean
        sigma2_upd = bcov[0, 0]
        Omega_upd = np.diag(bcov)[1:].copy()
        diff = np.max(np.abs(beta_upd - beta))
        beta0vec = bmat[:, 0].copy()
        beta1mat = bmat[:, 1:].copy()
        beta = beta_upd; sigma2 = sigma2_upd; OmegaDiag = Omega_upd
        if shrink:
            rho1, rho2 = _covlsconst_rho(scoreRaw, gg, Xf, Tvec, vcount,
                                         beta0vec, beta1mat, n, q)[:2]
            scoreWork = rho1 * gg + rho2 * scoreRaw
    return beta0vec, beta1mat, beta, sigma2, OmegaDiag, s


def _obj_func(score_flat, Xf, Tvec, vcount, beta0vec, beta1mat, beta0, beta1,
              sigma2, OmegaDiag, n, q):
    """Longitudinal objective (== obj_core / obj.func)."""
    start = np.concatenate([[0], np.cumsum(vcount)[:-1]])
    ll1 = 0.0
    for i in range(n):
        b_i = np.concatenate([[beta0vec[i]], beta1mat[i]])
        for v in range(vcount[i]):
            r = start[i] + v
            u1 = Xf[r] @ b_i
            ll1 += (u1 + score_flat[r] * np.exp(-u1)) * (Tvec[r] / 2.0)
    ll2 = np.sum(np.log(sigma2) / 2.0 + (beta0vec - beta0) ** 2 / (2.0 * sigma2))
    detO = np.prod(OmegaDiag)
    logdet = np.log(detO)
    Oginv = _ginv_diag(OmegaDiag)
    quad = ((beta1mat - beta1) ** 2) @ Oginv
    ll3 = np.sum(logdet / 2.0 + quad / 2.0)
    return ll1 + ll2 + ll3


def _sign_norm(g):
    g = g / np.sqrt(np.sum(g ** 2))
    if g[np.argmax(np.abs(g))] < 0:
        g = -g
    return g


# --------------------------------------------------------------------------- #
# given gamma: estimate the random-effects model                              #
# --------------------------------------------------------------------------- #
def cap_beta(Y, X, gamma, cov_shrinkage=True, max_itr=1000, tol=1e-4,
             score_return=True):
    """Estimate the random-effects model for a fixed ``gamma`` (== R ``cap_beta``)."""
    n = len(Y)
    p = Y[0][0].shape[1]
    q = np.asarray(X[0]).shape[1]
    gamma = np.asarray(gamma, dtype=float).ravel()
    gg = float(gamma @ gamma)

    nT_list = [np.array([Y[i][v].shape[0] for v in range(len(Y[i]))], dtype=float)
               for i in range(n)]
    if min(t.min() for t in nT_list) - 5 < p:
        cov_shrinkage = True

    if not cov_shrinkage:
        S, nT = _visit_covs(Y, shrinkage=False)
        sc = _scores_list(S, gamma)
        Xf, Tvec, vcount = _flatten(X, nT)
        scoreRaw = np.concatenate(sc)
        b0v, b1m, beta, sig2, Om = _init_re(sc, X, q)
        b0v, b1m, beta, sig2, Om, it = _re_recursion(
            scoreRaw, scoreRaw, gg, Xf, Tvec, vcount, False, n, q,
            b0v, b1m, beta, sig2, Om, max_itr, tol)
        score_out = sc
        shrink_par = None
    else:
        S_ls, nT = _visit_covs(Y, shrinkage=True)
        S_raw, _ = _visit_covs(Y, shrinkage=False)
        sc_ls = _scores_list(S_ls, gamma)
        sc_raw = _scores_list(S_raw, gamma)
        Xf, Tvec, vcount = _flatten(X, nT)
        scoreWork = np.concatenate(sc_ls)
        scoreRaw = np.concatenate(sc_raw)
        b0v, b1m, beta, sig2, Om = _init_re(sc_ls, X, q)
        b0v, b1m, beta, sig2, Om, it = _re_recursion(
            scoreWork, scoreRaw, gg, Xf, Tvec, vcount, True, n, q,
            b0v, b1m, beta, sig2, Om, max_itr, tol)
        rho1, rho2, mu, dlt, psi, phi = _covlsconst_rho(
            scoreRaw, gg, Xf, Tvec, vcount, b0v, b1m, n, q)
        scoreWork = rho1 * gg + rho2 * scoreRaw
        st = np.concatenate([[0], np.cumsum(vcount)[:-1]])
        score_out = [scoreWork[st[i]:st[i] + vcount[i]] for i in range(n)]
        shrink_par = {"rho1": rho1, "rho2": rho2, "mu": mu,
                      "phi2": phi, "psi2": psi, "delta2": dlt}

    res = {
        "gamma": gamma, "beta": beta,
        "beta0_random": b0v, "beta1_random": b1m,
        "beta0_sigma2": sig2,
        "beta1_Omega": Om,          # diagonal of the random-slope covariance
        "convergence": it < max_itr,
    }
    if shrink_par is not None:
        res["shrinkage"] = shrink_par
    if score_return:
        res["score"] = score_out
    return res


# --------------------------------------------------------------------------- #
# one direction (estimate gamma + beta) for a single gamma0                    #
# --------------------------------------------------------------------------- #
def _cap_d1_noshrink(Y, X, gamma0, max_itr, tol):
    """Joint gamma + random-effects update, no shrinkage (== capd1_core)."""
    n = len(Y)
    p = Y[0][0].shape[1]
    q = np.asarray(X[0]).shape[1]
    S, nT = _visit_covs(Y, shrinkage=False)
    Xf, Tvec, vcount = _flatten(X, nT)
    M = int(vcount.sum())
    Sig = np.empty((p, p, M))
    k = 0
    for i in range(n):
        for v in range(vcount[i]):
            Sig[:, :, k] = S[i][:, :, v]; k += 1
    H = (Sig * (Tvec / Tvec.sum())).sum(axis=2)

    gamma = np.asarray(gamma0, dtype=float).ravel()
    sc0 = _scores_list(S, gamma)
    b0v, b1m, beta, sig2, Om = _init_re(sc0, X, q)
    start = np.concatenate([[0], np.cumsum(vcount)[:-1]])

    s = 0
    diff = 100.0
    while s <= max_itr and diff > tol:
        s += 1
        score = np.einsum("i,ijr,j->r", gamma, Sig, gamma)
        dcov = np.concatenate([[sig2], Om])
        gdiag = _ginv_diag(dcov)
        bmat = np.empty((n, q))
        for i in range(n):
            b_i = np.concatenate([[b0v[i]], b1m[i]])
            pt1 = np.zeros(q); pt2 = np.zeros((q, q))
            for v in range(vcount[i]):
                r = start[i] + v
                xr = Xf[r]; u1 = xr @ b_i; e = np.exp(-u1); w = Tvec[r] / 2.0
                pt1 += (xr - score[r] * e * xr) * w
                pt2 += (score[r] * e) * np.outer(xr, xr) * w
            pt1 += gdiag * (b_i - beta)
            pt2[np.diag_indices(q)] += gdiag
            bmat[i] = b_i - np.linalg.solve(pt2, pt1)
        bmean = bmat.mean(axis=0)
        Bc = bmat - bmean
        bcov = (Bc.T @ Bc) / n
        beta_upd = bmean
        sig2_upd = bcov[0, 0]
        Om_upd = np.diag(bcov)[1:].copy()
        Amat = np.zeros((p, p))
        for i in range(n):
            bu = bmat[i]
            for v in range(vcount[i]):
                r = start[i] + v
                u1 = Xf[r] @ bu
                Amat += np.exp(-u1) * Sig[:, :, r] * (Tvec[r] / 2.0)
        gamma_upd = gamma_solve(Amat, H)
        diff = np.max(np.abs(beta_upd - beta))
        b0v = bmat[:, 0].copy(); b1m = bmat[:, 1:].copy()
        beta = beta_upd; sig2 = sig2_upd; Om = Om_upd
        gamma = gamma_upd
    return _sign_norm(gamma), s < max_itr


def _cap_d1_shrink(Y, X, gamma0, max_itr, tol):
    """Joint gamma + random-effects update, with shrinkage (== cap_D1 shrink path)."""
    n = len(Y)
    p = Y[0][0].shape[1]
    q = np.asarray(X[0]).shape[1]
    gamma = np.asarray(gamma0, dtype=float).ravel()
    gg = float(gamma @ gamma)

    S_ls, nT = _visit_covs(Y, shrinkage=True)
    S_raw, _ = _visit_covs(Y, shrinkage=False)
    Xf, Tvec, vcount = _flatten(X, nT)
    scoreRaw = np.concatenate(_scores_list(S_raw, gamma))
    start = np.concatenate([[0], np.cumsum(vcount)[:-1]])
    M = int(vcount.sum())
    raw_cube = np.empty((p, p, M)); k = 0
    for i in range(n):
        for v in range(vcount[i]):
            raw_cube[:, :, k] = S_raw[i][:, :, v]; k += 1
    eyeM = np.repeat(np.eye(p)[:, :, None], M, axis=2)

    sc_ls = _scores_list(S_ls, gamma)
    b0v, b1m, beta, sig2, Om = _init_re(sc_ls, X, q)
    score = np.concatenate(sc_ls)

    s = 0
    diff = 100.0
    while s <= max_itr and diff > tol:
        s += 1
        dcov = np.concatenate([[sig2], Om])
        gdiag = _ginv_diag(dcov)
        bmat = np.empty((n, q))
        for i in range(n):
            b_i = np.concatenate([[b0v[i]], b1m[i]])
            pt1 = np.zeros(q); pt2 = np.zeros((q, q))
            for v in range(vcount[i]):
                r = start[i] + v
                xr = Xf[r]; u1 = xr @ b_i; e = np.exp(-u1); w = Tvec[r] / 2.0
                pt1 += (xr - score[r] * e * xr) * w
                pt2 += (score[r] * e) * np.outer(xr, xr) * w
            pt1 += gdiag * (b_i - beta)
            pt2[np.diag_indices(q)] += gdiag
            bmat[i] = b_i - np.linalg.solve(pt2, pt1)
        bmean = bmat.mean(axis=0)
        Bc = bmat - bmean
        bcov = (Bc.T @ Bc) / n
        beta_upd = bmean
        sig2_upd = bcov[0, 0]
        Om_upd = np.diag(bcov)[1:].copy()
        b0v_upd = bmat[:, 0].copy(); b1m_upd = bmat[:, 1:].copy()
        rho1, rho2 = _covlsconst_rho(scoreRaw, gg, Xf, Tvec, vcount,
                                     b0v_upd, b1m_upd, n, q)[:2]
        ls_cube = rho1 * eyeM + rho2 * raw_cube
        score = rho1 * gg + rho2 * scoreRaw
        H = (ls_cube * (Tvec / Tvec.sum())).sum(axis=2)
        Amat = np.zeros((p, p))
        for i in range(n):
            bu = np.concatenate([[b0v_upd[i]], b1m_upd[i]])
            for v in range(vcount[i]):
                r = start[i] + v
                u1 = Xf[r] @ bu
                Amat += np.exp(-u1) * ls_cube[:, :, r] * (Tvec[r] / 2.0)
        gamma_upd = gamma_solve(Amat, H)
        diff = np.max(np.abs(beta_upd - beta))
        b0v = b0v_upd; b1m = b1m_upd
        beta = beta_upd; sig2 = sig2_upd; Om = Om_upd
        gamma = gamma_upd
    return _sign_norm(gamma), s < max_itr


def _cap_d1(Y, X, gamma0, cov_shrinkage, max_itr, tol):
    """One direction for a single gamma0: returns the re-estimated beta result."""
    if cov_shrinkage:
        g, conv = _cap_d1_shrink(Y, X, gamma0, max_itr, tol)
    else:
        g, conv = _cap_d1_noshrink(Y, X, gamma0, max_itr, tol)
    res = cap_beta(Y, X, g, cov_shrinkage=cov_shrinkage, max_itr=max_itr, tol=tol,
                   score_return=True)
    res["gamma"] = g
    res["convergence"] = conv
    return res


def _gamma0_mat(p, ninitial, seed):
    rng = np.random.RandomState(seed)
    M = rng.normal(size=(p, p + 1 + 5))
    M = M / np.sqrt((M ** 2).sum(axis=0, keepdims=True))
    cols = np.sort(rng.choice(M.shape[1], size=min(ninitial, M.shape[1]), replace=False))
    return M[:, cols]


def _cap_d1_opt(Y, X, cov_shrinkage, max_itr, tol, ninitial, seed):
    """Multi-initialization direction estimation; pick the lowest objective."""
    n = len(Y)
    p = Y[0][0].shape[1]
    q = np.asarray(X[0]).shape[1]
    if ninitial is None:
        ninitial = min(p + 1 + 5, 10)
    S, nT = _visit_covs(Y, shrinkage=cov_shrinkage)
    Xf, Tvec, vcount = _flatten(X, nT)
    M = int(vcount.sum())
    cube = np.empty((p, p, M)); k = 0
    for i in range(n):
        for v in range(vcount[i]):
            cube[:, :, k] = S[i][:, :, v]; k += 1
    H = (cube * (Tvec / Tvec.sum())).sum(axis=2)
    Sraw = S if not cov_shrinkage else _visit_covs(Y, shrinkage=False)[0]

    g0 = _gamma0_mat(p, ninitial, seed)
    best = None
    best_obj = np.inf
    for kk in range(g0.shape[1]):
        try:
            rk = _cap_d1(Y, X, g0[:, kk], cov_shrinkage, max_itr, tol)
        except Exception:
            continue
        g = rk["gamma"]
        g_un = g / np.sqrt(g @ H @ g)
        try:
            bt = cap_beta(Y, X, g_un, cov_shrinkage=cov_shrinkage, max_itr=max_itr,
                          tol=tol, score_return=False)
            if cov_shrinkage:
                raw = np.concatenate(_scores_list(Sraw, g_un))
                rho1, rho2 = _covlsconst_rho(raw, float(g_un @ g_un), Xf, Tvec,
                                             vcount, bt["beta0_random"],
                                             bt["beta1_random"], n, q)[:2]
                sflat = rho1 * float(g_un @ g_un) + rho2 * raw
            else:
                sflat = np.concatenate(_scores_list(S, g_un))
            Om = bt["beta1_Omega"]
            Omd = np.diag(Om) if np.ndim(Om) == 2 else np.atleast_1d(Om)
            obj = _obj_func(sflat, Xf, Tvec, vcount, bt["beta0_random"],
                            bt["beta1_random"], bt["beta"][0], bt["beta"][1:],
                            bt["beta0_sigma2"], Omd, n, q)
        except Exception:
            obj = np.inf
        if obj < best_obj:
            best_obj = obj
            best = rk
    return best


def _deflate(Y, Phi0, cov_shrinkage, beta_list, nT):
    """Project out identified directions (== cap_Dk deflation)."""
    n = len(Y)
    p = Y[0][0].shape[1]
    p0 = Phi0.shape[1]
    P = Phi0 @ Phi0.T
    Yt = []
    for i in range(n):
        vis = []
        for v in range(len(Y[i])):
            Yiv = np.asarray(Y[i][v], dtype=float)
            Y2 = Yiv - Yiv @ P
            if cov_shrinkage:
                vis.append(Y2)
            else:
                U, d, Vt = np.linalg.svd(Y2, full_matrices=False)
                b0 = np.array([beta_list[j]["beta0_random"][i] for j in range(p0)])
                dnew = d.copy()
                dnew[p - p0:p] = np.sqrt(np.exp(b0) * nT[i][v])
                vis.append((U * dnew) @ Vt)
        Yt.append(vis)
    return Yt


# --------------------------------------------------------------------------- #
# deviation-from-diagonality                                                   #
# --------------------------------------------------------------------------- #
def diag_level(Y, Phi, cov_shrinkage=True):
    """Per-direction deviation-from-diagonality (== R ``diag.level``)."""
    Phi = np.asarray(Phi, dtype=float)
    if Phi.ndim == 1 or Phi.shape[1] == 1:
        raise ValueError("Phi must have >= 2 columns")
    n = len(Y)
    ps = Phi.shape[1]
    S, nT = _visit_covs(Y, shrinkage=cov_shrinkage)
    Tsum = sum(t.sum() for t in nT)
    avg = np.ones(ps)
    for j in range(1, ps):
        prod = 1.0
        Pj = Phi[:, : j + 1]
        for i in range(n):
            for v in range(len(Y[i])):
                Mv = Pj.T @ S[i][:, :, v] @ Pj
                val = np.prod(np.diag(Mv)) / np.linalg.det(Mv)
                prod *= val ** (nT[i][v] / Tsum)
        avg[j] = prod
    return {"avg_level": avg}


# --------------------------------------------------------------------------- #
# main entry: find leading directions                                          #
# --------------------------------------------------------------------------- #
def cap_reg(Y, X, stop_crt="nD", nD=None, DfD_thred=5.0, OC=False,
            cov_shrinkage=True, max_itr=1000, tol=1e-4, ninitial=None,
            seed=500, verbose=False):
    """Longitudinal gamma-invariant CAP regression (== R ``capReg``).

    Note: the multi-start initial gamma matrix is drawn with numpy's RNG, so the
    specific starts differ from R's; the recovered optimum (direction + effects)
    matches because the objective minimum is start-robust.
    """
    n = len(Y)
    p = Y[0][0].shape[1]
    q = np.asarray(X[0]).shape[1]
    if stop_crt == "nD" and nD is None:
        stop_crt = "DfD"
    nT_list = [np.array([Y[i][v].shape[0] for v in range(len(Y[i]))], dtype=float)
               for i in range(n)]
    if min(t.min() for t in nT_list) - 5 < p:
        cov_shrinkage = True

    re1 = _cap_d1_opt(Y, X, cov_shrinkage, max_itr, tol, ninitial, seed)
    Phi = [re1["gamma"]]
    beta = [re1["beta"]]
    betas = [re1]
    if verbose:
        print("Component 1")

    def _grow(re_tmp):
        Phi.append(re_tmp["gamma"]); beta.append(re_tmp["beta"]); betas.append(re_tmp)
        if verbose:
            print(f"Component {len(Phi)}")

    if stop_crt == "nD":
        for j in range(2, nD + 1):
            Phi0 = np.column_stack(Phi)
            beta_list = [cap_beta(Y, X, Phi0[:, jj], cov_shrinkage=cov_shrinkage,
                                  max_itr=max_itr, tol=tol) for jj in range(Phi0.shape[1])]
            Yt = _deflate(Y, Phi0, cov_shrinkage, beta_list, nT_list)
            re_tmp = _cap_d1_opt(Yt, X, cov_shrinkage, max_itr, tol, ninitial, seed)
            if re_tmp is None:
                break
            _grow(re_tmp)
    elif stop_crt == "DfD":
        DfD = 1.0
        while DfD < DfD_thred:
            Phi0 = np.column_stack(Phi)
            beta_list = [cap_beta(Y, X, Phi0[:, jj], cov_shrinkage=cov_shrinkage,
                                  max_itr=max_itr, tol=tol) for jj in range(Phi0.shape[1])]
            Yt = _deflate(Y, Phi0, cov_shrinkage, beta_list, nT_list)
            re_tmp = _cap_d1_opt(Yt, X, cov_shrinkage, max_itr, tol, ninitial, seed)
            if re_tmp is None:
                break
            cand = np.column_stack(Phi + [re_tmp["gamma"]])
            DfD = diag_level(Y, cand, cov_shrinkage=cov_shrinkage)["avg_level"][-1]
            if DfD < DfD_thred and np.isfinite(DfD):
                _grow(re_tmp)
            else:
                break
    else:
        raise ValueError("stop_crt must be 'nD' or 'DfD'")

    Phi = np.column_stack(Phi)
    out = {
        "gamma": Phi,
        "beta": np.column_stack(beta),
        "beta0_random": np.column_stack([b["beta0_random"] for b in betas]),
        "beta0_sigma2": np.array([b["beta0_sigma2"] for b in betas]),
        "beta1_random": [b["beta1_random"] for b in betas],
        "beta1_Omega": [b["beta1_Omega"] for b in betas],
        "columns": [f"D{i + 1}" for i in range(Phi.shape[1])],
    }
    if Phi.shape[1] >= 2:
        out["DfD"] = diag_level(Y, Phi, cov_shrinkage=cov_shrinkage)
    return out


# --------------------------------------------------------------------------- #
# bootstrap inference                                                          #
# --------------------------------------------------------------------------- #
def cap_beta_boot(Y, X, gamma, cov_shrinkage=True, sims=1000,
                  boot_ci_type="perc", conf_level=0.95, boot_seed=100,
                  max_itr=1000, tol=1e-4, verbose=False):
    """Bootstrap inference for the fixed effects ``beta`` (== R ``cap_beta_boot``).

    Note: subjects are resampled with numpy's RNG (R uses ``set.seed`` per draw),
    so the specific resamples differ but the inference is statistically equivalent.
    """
    from scipy.stats import norm
    n = len(Y)
    gamma = np.asarray(gamma, dtype=float).ravel()
    base = cap_beta(Y, X, gamma, cov_shrinkage=cov_shrinkage, max_itr=max_itr, tol=tol)
    q = base["beta"].shape[0]
    boot = np.full((q, sims), np.nan)
    for b in range(sims):
        rng = np.random.RandomState(boot_seed + b)
        idx = rng.randint(0, n, size=n)
        Yt = [Y[i] for i in idx]
        Xt = [X[i] for i in idx]
        try:
            rt = cap_beta(Yt, Xt, gamma, cov_shrinkage=cov_shrinkage,
                          max_itr=max_itr, tol=tol)
            boot[:, b] = rt["beta"]
        except Exception:
            pass
        if verbose:
            print(f"Bootstrap sample {b + 1}")

    est = base["beta"]
    se = np.nanstd(boot, axis=1, ddof=1)
    stat = est / se
    pval = (1 - norm.cdf(np.abs(stat))) * 2
    lo, hi = (1 - conf_level) / 2, 1 - (1 - conf_level) / 2
    if boot_ci_type == "boot.se":
        z = norm.ppf(hi)
        ci = np.column_stack([est - z * se, est + z * se])
    else:                       # percentile
        ci = np.column_stack([np.nanquantile(boot, lo, axis=1),
                              np.nanquantile(boot, hi, axis=1)])
    return {
        "estimate": est, "se": se, "statistic": stat, "pvalue": pval,
        "lower": ci[:, 0], "upper": ci[:, 1], "boot": boot,
    }
