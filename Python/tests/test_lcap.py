"""Self-contained tests for the LCAP (gamma-invariant longitudinal) port."""
import numpy as np
from capcov import lcap


def _simulate(n=40, p=6, nV=6, Ti=60, seed=2024):
    rng = np.random.RandomState(seed)
    tg = np.linspace(-0.5, 0.5, nV)
    A = rng.normal(size=(p, p))
    Phi, _ = np.linalg.qr(A)
    if Phi[np.argmax(np.abs(Phi[:, 0])), 0] < 0:
        Phi[:, 0] = -Phi[:, 0]
    base = np.full(p, 0.5)
    Y, X = [], []
    truth_b1 = np.array([-0.5, 0.6])     # (time, dose) fixed slopes
    for i in range(n):
        b0 = rng.normal(0, 0.3)
        bt = truth_b1[0] + rng.normal(0, 0.2)
        bd = truth_b1[1] + rng.normal(0, 0.2)
        dose = np.round(rng.normal(size=nV), 2)
        Xi = np.column_stack([np.ones(nV), tg, dose])
        vis = []
        for v in range(nV):
            e = base.copy(); e[0] = np.exp(b0 + bt * tg[v] + bd * dose[v])
            S = Phi @ np.diag(e) @ Phi.T; S = (S + S.T) / 2
            w, V = np.linalg.eigh(S)
            rt = V @ np.diag(np.sqrt(np.clip(w, 1e-8, None))) @ V.T
            vis.append(rng.normal(size=(Ti, p)) @ rt)
        Y.append(vis); X.append(Xi)
    return dict(Y=Y, X=X, gamma=Phi[:, 0], beta1=truth_b1)


def test_cap_beta_recovers_fixed_effects_no_shrinkage():
    d = _simulate()
    r = lcap.cap_beta(d["Y"], d["X"], d["gamma"], cov_shrinkage=False)
    assert r["beta"].shape[0] == 3
    assert abs(r["beta"][1] - d["beta1"][0]) < 0.1     # time
    assert abs(r["beta"][2] - d["beta1"][1]) < 0.1     # dose
    assert r["beta0_sigma2"] > 0
    assert np.all(np.atleast_1d(r["beta1_Omega"]) >= 0)


def test_cap_beta_shrinkage_runs():
    d = _simulate()
    r = lcap.cap_beta(d["Y"], d["X"], d["gamma"], cov_shrinkage=True)
    assert "shrinkage" in r
    assert np.isfinite(r["beta"]).all()
    assert abs(r["beta"][1] - d["beta1"][0]) < 0.15


def test_cap_d1_recovers_direction_no_shrinkage():
    d = _simulate()
    rng = np.random.RandomState(7)
    g0 = rng.normal(size=d["gamma"].shape[0]); g0 /= np.linalg.norm(g0)
    r = lcap._cap_d1(d["Y"], d["X"], g0, False, 1000, 1e-4)
    cos = abs(r["gamma"] @ d["gamma"])
    assert cos > 0.99


def test_cap_reg_direction_and_shapes():
    d = _simulate()
    fit = lcap.cap_reg(d["Y"], d["X"], stop_crt="nD", nD=1, cov_shrinkage=False,
                       ninitial=5)
    assert fit["gamma"].shape == (6, 1)
    assert abs(fit["gamma"][:, 0] @ d["gamma"]) > 0.99
    assert fit["beta"].shape[0] == 3


def test_diag_level_two_directions():
    d = _simulate()
    Phi = np.linalg.qr(np.random.RandomState(3).normal(size=(6, 6)))[0][:, :2]
    dl = lcap.diag_level(d["Y"], Phi, cov_shrinkage=False)
    assert dl["avg_level"].shape[0] == 2
    assert dl["avg_level"][0] == 1.0
    assert dl["avg_level"][1] >= 1.0          # deviation-from-diagonality >= 1


def test_cap_beta_boot_inference():
    d = _simulate()
    bi = lcap.cap_beta_boot(d["Y"], d["X"], d["gamma"], cov_shrinkage=False, sims=30)
    assert bi["estimate"].shape[0] == 3
    assert np.all(np.isfinite(bi["se"]))
