"""Self-contained tests for the CAP-CoC port (no R dependency)."""
import numpy as np
import capcov
from capcov import coc


def _simulate(n=50, p=12, q=10, Tx=80, Ty=80, seed=1):
    rng = np.random.RandomState(seed)
    W = np.column_stack([np.ones(n), np.tile([0.0, 1.0], n // 2 + 1)[:n], rng.normal(size=n)])
    W[:, 2] = (W[:, 2] - W[:, 2].mean()) / W[:, 2].std(ddof=1)
    alpha, beta = 0.8, np.array([0.2, 0.5, -0.3])
    Phix, _ = np.linalg.qr(rng.normal(size=(p, p)))
    Phiy, _ = np.linalg.qr(rng.normal(size=(q, q)))
    X, Y = [], []
    for i in range(n):
        lx = rng.normal()
        evx = np.full(p, 0.5); evx[0] = np.exp(lx)
        Sx = Phix @ np.diag(evx) @ Phix.T; Sx = (Sx + Sx.T) / 2
        w, V = np.linalg.eigh(Sx); rtx = V @ np.diag(np.sqrt(np.clip(w, 1e-8, None))) @ V.T
        X.append(rng.normal(size=(Tx, p)) @ rtx)
        ly = alpha * lx + beta @ W[i]
        evy = np.full(q, 0.4); evy[0] = np.exp(ly)
        Sy = Phiy @ np.diag(evy) @ Phiy.T; Sy = (Sy + Sy.T) / 2
        w, V = np.linalg.eigh(Sy); rty = V @ np.diag(np.sqrt(np.clip(w, 1e-8, None))) @ V.T
        Y.append(rng.normal(size=(Ty, q)) @ rty)
    return dict(Y=Y, X=X, W=W, gamma=Phiy[:, 0], theta=Phix[:, 0], alpha=alpha, beta=beta)


def test_cov_sk_cubes_symmetric_pd():
    d = _simulate()
    skx = coc._cov_sk_x(d["X"])
    sky = coc._cov_sk_y(d["Y"], d["gamma"], np.ones(len(d["Y"])))
    for cube in (skx, sky):
        S0 = cube[:, :, 0]
        assert np.allclose(S0, S0.T)
        assert np.all(np.linalg.eigvalsh(S0) > 0)     # shrinkage -> PD


def test_coef_recovers_truth_no_shrinkage():
    d = _simulate()
    c = coc.coc_coef(d["Y"], d["X"], d["W"], d["gamma"], d["theta"],
                     cov_shrinkage_y=False, cov_shrinkage_x=False, max_itr=1000)
    assert abs(c["alpha"] - d["alpha"]) < 0.15
    assert np.max(np.abs(c["beta"] - d["beta"])) < 0.15


def test_coef_equals_joint_ols_no_shrinkage():
    # with fixed scores the alternating update converges to the joint OLS fit
    d = _simulate()
    Sx = coc._sample_cov_cube(d["X"]); Sy = coc._sample_cov_cube(d["Y"])
    lsx = np.log(coc.scores(Sx, d["theta"]))
    lsy = np.log(coc.scores(Sy, d["gamma"]))
    D = np.column_stack([lsx, d["W"]])
    coef_ols = np.linalg.lstsq(D, lsy, rcond=None)[0]
    c = coc.coc_coef(d["Y"], d["X"], d["W"], d["gamma"], d["theta"],
                     cov_shrinkage_y=False, cov_shrinkage_x=False, max_itr=2000, tol=1e-8)
    assert abs(c["alpha"] - coef_ols[0]) < 1e-4
    assert np.max(np.abs(c["beta"] - coef_ols[1:])) < 1e-4


def test_coc_reg_two_directions_shapes_and_orthogonality():
    d = _simulate()
    fit = coc.coc_reg(d["Y"], d["X"], d["W"], stop_crt="nD", nD=2,
                      burn_in=50, max_itr=300, ninitial=5)
    assert fit["gamma"].shape == (10, 2)
    assert fit["theta"].shape == (12, 2)
    assert fit["beta"].shape == (3, 2)
    assert fit["alpha"].shape[0] == 2
    assert fit["DfD_y"]["avg_level"].shape[0] == 2
    # each gamma/theta column is unit norm
    assert np.allclose(np.linalg.norm(fit["gamma"], axis=0), 1.0)
    assert np.allclose(np.linalg.norm(fit["theta"], axis=0), 1.0)


def test_asymptotic_inference_finite():
    d = _simulate()
    inf = coc.coc_coef_asmp(d["Y"], d["X"], d["W"], d["gamma"], d["theta"],
                            cov_shrinkage_y=False, cov_shrinkage_x=False)
    assert np.isfinite(inf["alpha"]["se"])
    assert np.all(np.isfinite(inf["beta"]["se"]))
    assert inf["beta"]["estimate"].shape[0] == 3


def test_bootstrap_inference_runs():
    d = _simulate()
    bi = coc.coc_coef_boot(d["Y"], d["X"], d["W"], d["gamma"], d["theta"],
                           cov_shrinkage_y=False, cov_shrinkage_x=False, sims=40)
    assert np.isfinite(bi["alpha"]["se"])
    assert bi["beta"]["estimate"].shape[0] == 3
