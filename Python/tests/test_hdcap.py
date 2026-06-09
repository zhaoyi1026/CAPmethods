"""Self-contained tests for the HDCAP / CAP port (no R dependency)."""
import numpy as np
import capcov
from capcov import _core


def _simulate(n=60, p=8, Ti=120, seed=11):
    rng = np.random.RandomState(seed)
    group = np.tile([0.0, 1.0], n // 2 + 1)[:n]
    age = rng.normal(size=n); age = (age - age.mean()) / age.std(ddof=1)
    X = np.column_stack([np.ones(n), group, age])
    A = rng.normal(size=(p, p)); Phi, _ = np.linalg.qr(A)
    if Phi[np.argmax(np.abs(Phi[:, 0])), 0] < 0:
        Phi[:, 0] = -Phi[:, 0]
    b1 = np.array([0.0, 0.9, -0.5])
    base = np.array([np.nan, .6, .45, .3, .25, .2, .15, .1])[:p]
    Y = []
    for i in range(n):
        ev = base.copy(); ev[0] = np.exp(X[i] @ b1)
        S = Phi @ np.diag(ev) @ Phi.T; S = (S + S.T) / 2
        w, V = np.linalg.eigh(S)
        root = V @ np.diag(np.sqrt(np.clip(w, 1e-8, None))) @ V.T
        Y.append(rng.normal(size=(Ti, p)) @ root)
    return Y, X, Phi[:, 0], b1


def test_cov_ls_symmetric():
    rng = np.random.RandomState(0)
    Z = rng.normal(size=(40, 6))
    S = _core.cov_ls(Z)
    assert np.allclose(S, S.T)
    assert np.all(np.linalg.eigvalsh(S) > 0)        # shrinkage -> PD


def test_scores_and_accum():
    rng = np.random.RandomState(1)
    Y, X, _, _ = _simulate()
    Sigma, Tvec = _core.sigma_cube(Y)
    v = rng.normal(size=Sigma.shape[0])
    ref = np.array([v @ Sigma[:, :, i] @ v for i in range(Sigma.shape[2])])
    assert np.allclose(_core.scores(Sigma, v), ref)
    w = rng.uniform(size=Sigma.shape[2])
    ref_acc = sum(w[i] * Sigma[:, :, i] for i in range(Sigma.shape[2]))
    assert np.allclose(_core.accum(Sigma, w), ref_acc)


def test_cap_recovers_truth_no_shrinkage():
    Y, X, gamma_true, beta_true = _simulate()
    res = capcov.cap_reg(Y, X, stop_crt="nD", nD=1, cov_shrinkage=False)
    g = res["gamma"][:, 0]
    cos = abs(g @ gamma_true) / (np.linalg.norm(g) * np.linalg.norm(gamma_true))
    assert cos > 0.99
    assert abs(res["beta"][1, 0] - beta_true[1]) < 0.1   # group
    assert abs(res["beta"][2, 0] - beta_true[2]) < 0.1   # age


def test_hdcap_shrinkage_runs_and_recovers():
    Y, X, gamma_true, beta_true = _simulate()
    res = capcov.cap_reg(Y, X, stop_crt="nD", nD=1, cov_shrinkage=True)
    g = res["gamma"][:, 0]
    cos = abs(g @ gamma_true) / (np.linalg.norm(g) * np.linalg.norm(gamma_true))
    assert cos > 0.99
    assert "shrinkage" in res


def test_dfd_two_directions():
    Y, X, _, _ = _simulate()
    res = capcov.cap_reg(Y, X, stop_crt="nD", nD=2, cov_shrinkage=False)
    assert res["gamma"].shape[1] == 2
    assert res["DfD"]["avg_level"].shape[0] == 2


def test_bootstrap_inference():
    Y, X, gamma_true, _ = _simulate()
    res = capcov.cap_reg(Y, X, stop_crt="nD", nD=1, cov_shrinkage=False)
    bi = capcov.cap_beta_boot(Y, X, res["gamma"][:, 0], cov_shrinkage=False, sims=50)
    assert bi["estimate"].shape[0] == X.shape[1]
    assert np.all(np.isfinite(bi["se"]))
