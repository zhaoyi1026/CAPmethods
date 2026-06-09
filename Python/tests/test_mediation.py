"""Self-contained tests for the CAP-mediation port."""
import numpy as np
from capcov import mediation as md


def _simulate(n=60, p=12, Ti=60, seed=1):
    rng = np.random.RandomState(seed)
    x = np.tile([0.0, 1.0], n // 2 + 1)[:n]
    w = rng.normal(size=n); w = (w - w.mean()) / w.std(ddof=1)
    X = np.column_stack([x, w])
    a0, ax, aw = 0.2, 0.8, -0.3
    beta, g0, gx, gw = 0.7, 0.1, 0.4, 0.2
    Phi = np.linalg.qr(rng.normal(size=(p, p)))[0]
    if Phi[np.argmax(np.abs(Phi[:, 0])), 0] < 0:
        Phi[:, 0] = -Phi[:, 0]
    base = np.concatenate([[np.nan], np.exp(np.linspace(np.log(1.2), np.log(0.3), p - 1))])
    M, Y = [], np.empty(n)
    for i in range(n):
        lv = a0 + ax * x[i] + aw * w[i] + rng.normal(0, 0.3)
        ev = base.copy(); ev[0] = np.exp(lv)
        S = Phi @ np.diag(ev) @ Phi.T; S = (S + S.T) / 2
        wv, V = np.linalg.eigh(S)
        rt = V @ np.diag(np.sqrt(np.clip(wv, 1e-8, None))) @ V.T
        M.append(rng.normal(size=(Ti, p)) @ rt)
        Y[i] = g0 + gx * x[i] + gw * w[i] + beta * lv + rng.normal(0, 0.3)
    return dict(X=X, M=M, Y=Y, theta=Phi[:, 0], IE=ax * beta)


def test_med_cov_symmetric():
    d = _simulate(n=10, Ti=40)
    S, nT = md._med_cov(d["M"])
    assert np.allclose(S[:, :, 0], S[:, :, 0].T)
    assert nT[0] == 40


def test_coef_identifiable_outputs_robust_to_blup():
    d = _simulate()
    c0 = md.cap_med_coef(d["X"], d["M"], d["Y"], d["theta"], blup_shrink=0.0)
    c1 = md.cap_med_coef(d["X"], d["M"], d["Y"], d["theta"], blup_shrink=1.0)
    # alpha / beta / gamma / IE are insensitive to the (non-identifiable) BLUP split
    assert np.allclose(c0["alpha"], c1["alpha"])
    assert abs(c0["beta"] - c1["beta"]) < 1e-10
    assert abs(c0["IE"] - c1["IE"]) < 1e-10


def test_d1_recovers_direction_and_effects():
    d = _simulate()
    f = md.cap_med_d1(d["X"], d["M"], d["Y"], max_itr=500)
    assert abs(f["theta"] @ d["theta"]) > 0.95           # recovers projection
    assert abs(f["alpha"][0] - 0.8) < 0.15               # exposure -> mediator
    assert f["IE"] > 0.2                                  # positive indirect effect


def test_cap_med_main_shapes():
    d = _simulate()
    fit = md.cap_med(d["X"], d["M"], d["Y"], stop_crt="nD", nD=1, ninitial=3)
    assert fit["theta"].shape == (12, 1)
    assert fit["alpha"].shape == (2, 1)
    assert fit["beta"].shape[0] == 1
    assert abs(fit["theta"][:, 0] @ d["theta"]) > 0.95
