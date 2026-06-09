"""Self-contained tests for the CAP-clustering (PCL) port."""
import numpy as np
from capcov import clustering as cl


def _simulate(n=120, p=5, Ti=200, seed=7):
    rng = np.random.RandomState(seed)
    W = np.column_stack([np.ones(n), rng.normal(size=n)])
    clt = (1.0 / (1.0 + np.exp(-(W @ np.array([0.0, 3.0])))) > rng.uniform(size=n)).astype(int)
    X = np.column_stack([np.ones(n), np.tile([0.0, 1.0], n // 2 + 1)[:n]])
    beta_true = np.array([[0.0, 1.5], [2.0, -1.5]])    # q1=2 x K=2, well separated
    Phi = np.linalg.qr(rng.normal(size=(p, p)))[0]
    base = np.linspace(0.8, 0.1, p)
    Y = []
    for i in range(n):
        ev = base.copy(); ev[0] = np.exp(X[i] @ beta_true[:, clt[i]])
        S = Phi @ np.diag(ev) @ Phi.T; S = (S + S.T) / 2
        w, V = np.linalg.eigh(S)
        Y.append(rng.normal(size=(Ti, p)) @ (V @ np.diag(np.sqrt(np.clip(w, 1e-8, None))) @ V.T))
    return dict(Y=Y, X=X, W=W, gamma=Phi[:, 0], cl=clt, beta=beta_true)


def test_kernels_second_moment_and_score():
    d = _simulate(n=20, Ti=50)
    S, Tvec = cl._smat(d["Y"])
    g = d["gamma"]
    Yi = d["Y"][3]
    assert abs(cl._score(S, g)[3] - ((Yi @ g) ** 2).sum() / Yi.shape[0]) < 1e-12
    assert np.allclose(S[:, :, 0], S[:, :, 0].T)


def test_coef_recovers_truth_well_separated():
    d = _simulate()
    f = cl.cap_pcl_coef(d["Y"], d["X"], d["W"], d["gamma"], ncluster=2, nstart=10)
    agree = max(np.mean((f["class"] - 1) == d["cl"]),
                np.mean((f["class"] - 1) == (1 - d["cl"])))
    assert agree > 0.95
    perm = cl._align_clusters(d["beta"], f["beta"])
    assert np.max(np.abs(f["beta"][:, perm] - d["beta"])) < 0.2


def test_coef_shapes_and_loglik():
    d = _simulate(n=40, Ti=80)
    f = cl.cap_pcl_coef(d["Y"], d["X"], d["W"], d["gamma"], ncluster=2, nstart=3)
    assert f["beta"].shape == (2, 2)
    assert f["alpha"].shape == (2, 2)
    assert f["class"].shape[0] == 40
    assert np.isfinite(f["logLik"])


def test_diag_level_matches_manual():
    d = _simulate(n=15, Ti=60)
    G = np.linalg.qr(np.random.RandomState(0).normal(size=(5, 5)))[0][:, :2]
    dl = cl.diag_level(d["Y"], G)
    S, _ = cl._smat(d["Y"])
    M = G.T @ S[:, :, 0] @ G
    manual = np.prod(np.diag(M)) / np.linalg.det(M)
    assert abs(dl["sub_level"][0, 1] - manual) < 1e-12
    assert dl["avg_level"][0] == 1.0


def test_firth_logistic_runs():
    rng = np.random.RandomState(0)
    W = np.column_stack([np.ones(50), rng.normal(size=50)])
    y = (W @ np.array([0.0, 1.5]) + rng.normal(size=50) > 0).astype(float)
    coef = cl._firth_logistic_weighted(W, y, np.ones(50))
    assert coef.shape[0] == 2 and np.all(np.isfinite(coef))
