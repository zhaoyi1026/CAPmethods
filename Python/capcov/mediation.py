"""CAP-mediation -- Python port of CAP-mediation/V4.

Covariance-matrix *mediator* model: a scalar exposure ``X`` affects a covariance-
valued mediator ``M`` (summarized along a projection ``theta``) which affects a
scalar outcome ``Y``:

    M model:  log(theta' Sigma_i theta) = alpha0 + X_i' alpha + b_i      (random intercept b_i)
    Y model:  Y_i = gamma0 + X_i' gamma + beta * log(theta' Sigma_i theta) + e_i

The indirect (mediation) effect is ``IE = alpha_x * beta``. ``cap_med`` finds the
projection theta and the path coefficients.

Ported from ``CAP-mediation/V4/CAPMediation.R`` and ``med_kernels.cpp``.

Important fidelity note
-----------------------
R fits the M model with ``nlme::lme(score ~ X, random = ~1|ID)`` where each
subject contributes a **single** score, i.e. the random-intercept groups are
singletons. The random-intercept / residual variance split is then
*non-identifiable* (the marginal likelihood depends only on their sum), so R's
``lme`` returns an optimizer-arbitrary split and BLUP ``alpha0.rnd``. This port
uses a self-contained, deterministic random-intercept fit (fixed effects = GLS =
OLS, which are identifiable; ``blup_shrink`` controls the BLUP). The identifiable
quantities -- ``alpha`` (M-model), ``beta`` (Y-model), the projection ``theta``,
and ``IE`` -- agree with R; the non-identifiable random effect may differ.
statsmodels is intentionally avoided (not required; it segfaults in some
environments). See the package README.
"""
from __future__ import annotations

import numpy as np

from ._core import gamma_solve

__all__ = ["cap_med", "cap_med_coef", "cap_med_d1"]


# --------------------------------------------------------------------------- #
# kernels (== med_kernels.cpp)                                                #
# --------------------------------------------------------------------------- #
def _med_cov(M):
    """Per-subject sample covariance cube ``cov(M_i)`` (denominator T_i - 1)."""
    n = len(M)
    p = np.asarray(M[0]).shape[1]
    S = np.empty((p, p, n))
    nT = np.empty(n)
    for i in range(n):
        Mi = np.asarray(M[i], dtype=float)
        T = Mi.shape[0]
        nT[i] = T
        Mc = Mi - Mi.mean(axis=0, keepdims=True)
        S[:, :, i] = (Mc.T @ Mc) / (T - 1)
    return S, nT


def _score(S, v):
    v = np.asarray(v, dtype=float).ravel()
    return np.einsum("i,ijk,j->k", v, S, v)


def _accum(S, w):
    return S @ np.asarray(w, dtype=float).ravel()


# --------------------------------------------------------------------------- #
# given theta: estimate path coefficients (== CAPMediation_coef)               #
# --------------------------------------------------------------------------- #
def cap_med_coef(X, M, Y, theta, M_cov=None, nT=None, blup_shrink=1.0):
    """Estimate the mediation path coefficients for a fixed ``theta``.

    ``blup_shrink`` in [0, 1] sets the random-intercept BLUP weight (1 = each
    subject's own log-variance, the natural CAP/saturated choice; 0 = no random
    effect). The identifiable outputs (alpha, beta, gamma, theta, IE) are
    insensitive to it; only ``alpha0_rnd`` / ``tau2`` depend on it.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    n = len(Y)
    if M_cov is None:
        M_cov, nT = _med_cov(M)
    theta = np.asarray(theta, dtype=float).ravel()
    score = _score(M_cov, theta)
    logs = np.log(score)

    # Y model: Y ~ [1, X, log(score)]  (OLS / pinv, == ginv solve)
    Z = np.column_stack([np.ones(n), X, logs])
    mu = np.linalg.pinv(Z.T @ Z) @ (Z.T @ Y)
    gamma0 = mu[0]
    gamma = mu[1:-1]
    beta = mu[-1]
    sigma2 = np.mean((Y - Z @ mu) ** 2)

    # M model fixed effects: GLS == OLS of log(score) on [1, X] (singleton groups)
    Zm = np.column_stack([np.ones(n), X])
    am = np.linalg.pinv(Zm.T @ Zm) @ (Zm.T @ logs)
    alpha0 = am[0]
    alpha = am[1:]
    fitted = Zm @ am
    # alpha0_rnd is the per-subject *intercept* (fixed intercept + shrunk random
    # intercept BLUP); the X-slope effect is added separately as X @ alpha in the
    # linear predictor, so it must NOT be included here (matches R's lme
    # coefficients$random$ID + fixed[1]).
    alpha0_rnd = alpha0 + blup_shrink * (logs - fitted)
    tau2 = np.mean((alpha0_rnd - alpha0) ** 2)

    IE = alpha[0] * beta
    return {
        "theta": theta, "alpha": alpha, "beta": float(beta), "gamma": gamma,
        "IE": float(IE), "alpha0": float(alpha0), "alpha0_rnd": alpha0_rnd,
        "gamma0": float(gamma0), "tau2": float(tau2), "sigma2": float(sigma2),
        "score": score,
    }


def _obj_func(X, M_cov, Y, nT, theta, c):
    """Mediation objective (== obj.func)."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    score = _score(M_cov, theta)
    a0r = c["alpha0_rnd"]; alpha = c["alpha"]
    lp = a0r + X @ alpha
    ll1 = np.sum((lp + score * np.exp(-lp)) * nT) / 2.0
    ll2 = np.sum((Y - c["gamma0"] - X @ c["gamma"] - c["beta"] * np.log(score)) ** 2
                 / c["sigma2"] + np.log(c["sigma2"])) / 2.0
    ll3 = np.sum((a0r - c["alpha0"]) ** 2 / c["tau2"] + np.log(c["tau2"])) / 2.0
    return ll1 + ll2 + ll3


# --------------------------------------------------------------------------- #
# one direction: estimate theta + coefficients (== CAPMediation_D1)            #
# --------------------------------------------------------------------------- #
def _sign_norm(v):
    v = v / np.sqrt(np.sum(v ** 2))
    if v[np.argmax(np.abs(v))] < 0:
        v = -v
    return v


def cap_med_d1(X, M, Y, H=None, theta0=None, max_itr=1000, tol=1e-4,
               blup_shrink=1.0, M_cov=None, nT=None):
    """Estimate theta and the path coefficients (== R ``CAPMediation_D1``)."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    if M_cov is None:
        M_cov, nT = _med_cov(M)
    p = M_cov.shape[0]
    if H is None:
        H = M_cov.mean(axis=2)
    if theta0 is None:
        theta0 = np.full(p, 1.0 / np.sqrt(p))
    theta0 = np.asarray(theta0, dtype=float).ravel()

    c = cap_med_coef(X, M, Y, theta0, M_cov=M_cov, nT=nT, blup_shrink=blup_shrink)
    s = 0
    diff = 100.0
    while s <= max_itr and diff > tol:
        s += 1
        score = _score(M_cov, theta0)
        U = np.exp(-c["alpha0_rnd"] - X @ c["alpha"])
        V = Y - c["gamma0"] - X @ c["gamma"]
        beta = c["beta"]
        wA = nT * U - (2 * beta * (V - beta * np.log(score))) / (c["sigma2"] * score)
        A = _accum(M_cov, wA)
        theta_new = gamma_solve(A, H)
        c_new = cap_med_coef(X, M, Y, theta_new, M_cov=M_cov, nT=nT, blup_shrink=blup_shrink)
        diff = np.max(np.abs(np.concatenate([
            c_new["alpha"] - c["alpha"], [c_new["beta"] - c["beta"]],
            c_new["gamma"] - c["gamma"]])))
        c = c_new
        theta0 = theta_new

    theta_est = _sign_norm(theta0)
    out = cap_med_coef(X, M, Y, theta_est, M_cov=M_cov, nT=nT, blup_shrink=blup_shrink)
    out["obj"] = _obj_func(X, M_cov, Y, nT, theta_est, out)
    out["convergence"] = s < max_itr
    return out


def _theta0_mat(p, ninitial, seed):
    rng = np.random.RandomState(seed)
    Tm = rng.normal(size=(p, p + 1 + 5))
    Tm = Tm / np.sqrt((Tm ** 2).sum(axis=0, keepdims=True))
    cols = np.sort(rng.choice(Tm.shape[1], size=min(ninitial, Tm.shape[1]), replace=False))
    return Tm[:, cols]


def _cap_med_d1_opt(X, M, Y, M_cov, nT, H, max_itr, tol, ninitial, seed, blup_shrink):
    p = M_cov.shape[0]
    if ninitial is None:
        ninitial = min(p, 10)
    Tmat = _theta0_mat(p, ninitial, seed)
    best = None
    best_obj = np.inf
    for kk in range(Tmat.shape[1]):
        try:
            rk = cap_med_d1(X, M, Y, H=H, theta0=Tmat[:, kk], max_itr=max_itr,
                            tol=tol, blup_shrink=blup_shrink, M_cov=M_cov, nT=nT)
        except Exception:
            continue
        # Select on the objective evaluated at the H-UNSCALED theta (matches R's
        # CAPMediation_D1_opt). The 1/sqrt(theta'H theta) rescaling penalizes
        # high-variance background directions, so the mediating direction wins.
        try:
            g = rk["theta"]
            g_un = g / np.sqrt(g @ H @ g)
            c_un = cap_med_coef(X, M, Y, g_un, M_cov=M_cov, nT=nT, blup_shrink=blup_shrink)
            obj = _obj_func(X, M_cov, Y, nT, g_un, c_un)
        except Exception:
            obj = np.inf
        if obj < best_obj:
            best_obj = obj
            best = rk
    return best


def _diag_level(M_cov, Theta, nT):
    Theta = np.asarray(Theta, dtype=float)
    if Theta.ndim == 1 or Theta.shape[1] == 1:
        return {"avg_level": np.array([1.0])}
    n = M_cov.shape[2]
    r = Theta.shape[1]
    dl = np.ones((n, r))
    for i in range(n):
        for j in range(1, r):
            P = Theta[:, : j + 1]
            Mt = P.T @ M_cov[:, :, i] @ P
            dl[i, j] = np.prod(np.diag(Mt)) / np.linalg.det(Mt)
    avg = np.array([np.prod(dl[:, j] ** (nT / nT.sum())) for j in range(r)])
    return {"avg_level": avg, "sub_level": dl}


# --------------------------------------------------------------------------- #
# main entry: leading mediation directions (== CAPMediation)                    #
# --------------------------------------------------------------------------- #
def cap_med(X, M, Y, H=None, stop_crt="nD", nD=None, DfD_thred=2.0, Y_remove=True,
            max_itr=1000, tol=1e-4, ninitial=None, seed=100, blup_shrink=1.0,
            verbose=False):
    """Covariance-mediator CAP regression (== R ``CAPMediation``).

    Returns ``theta`` (p x k), ``alpha`` (path a), ``beta`` (path b), ``gamma``
    (direct effect), and ``IE`` (indirect effect) per direction.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    M_cov, nT = _med_cov(M)
    if stop_crt == "nD" and nD is None:
        stop_crt = "DfD"
    if H is None:
        H = M_cov.mean(axis=2)

    re1 = _cap_med_d1_opt(X, M, Y, M_cov, nT, H, max_itr, tol, ninitial, seed, blup_shrink)
    Theta = [re1["theta"]]
    res = [re1]
    if verbose:
        print("Component 1")

    def _grow(rt):
        Theta.append(rt["theta"]); res.append(rt)
        if verbose:
            print(f"Component {len(Theta)}")

    def _deflate(Theta0):
        P = Theta0 @ Theta0.T
        Mt = [np.asarray(M[i], dtype=float) - np.asarray(M[i], dtype=float) @ P
              for i in range(len(M))]
        Sc, nt = _med_cov(Mt)
        return Mt, Sc, nt

    if stop_crt == "nD":
        for _ in range(2, nD + 1):
            Mt, Sc, nt = _deflate(np.column_stack(Theta))
            rt = _cap_med_d1_opt(X, Mt, Y, Sc, nt, Sc.mean(axis=2), max_itr, tol,
                                 ninitial, seed, blup_shrink)
            if rt is None:
                break
            _grow(rt)
    elif stop_crt == "DfD":
        DfD = 1.0
        while DfD < DfD_thred:
            Mt, Sc, nt = _deflate(np.column_stack(Theta))
            rt = _cap_med_d1_opt(X, Mt, Y, Sc, nt, Sc.mean(axis=2), max_itr, tol,
                                 ninitial, seed, blup_shrink)
            if rt is None:
                break
            cand = np.column_stack(Theta + [rt["theta"]])
            DfD = _diag_level(M_cov, cand, nT)["avg_level"][-1]
            if DfD < DfD_thred and np.isfinite(DfD):
                _grow(rt)
            else:
                break
    else:
        raise ValueError("stop_crt must be 'nD' or 'DfD'")

    Tmat = np.column_stack(Theta)
    out = {
        "theta": Tmat,
        "alpha": np.column_stack([r["alpha"] for r in res]),
        "beta": np.array([r["beta"] for r in res]),
        "gamma": np.column_stack([r["gamma"] for r in res]),
        "IE": np.array([r["IE"] for r in res]),
        "columns": [f"C{i + 1}" for i in range(Tmat.shape[1])],
    }
    if Tmat.shape[1] >= 2:
        out["DfD"] = _diag_level(M_cov, Tmat, nT)
    return out
