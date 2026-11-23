"""Tests of the 0.2.0 changes (as HDCAP V9 / CAP-CoC V4 of the R package CAPmethods 0.2.0).

* R's random numbers for the default initial value of beta;
* the objective that compares initial values is the function the fit minimizes;
* several initial values of beta (``beta0``), the converged fit with the lowest
  objective kept;
* a first direction that cannot be estimated is a clear error, a later one ends
  the search.
"""
import contextlib

import numpy as np

import capcov
from capcov import _core, hdcap, coc, lcap, mediation, clustering
from test_hdcap import _simulate as _sim_cap
from test_coc import _simulate as _sim_coc
from test_lcap import _simulate as _sim_lcap
from test_mediation import _simulate as _sim_med
from test_clustering import _simulate as _sim_pcl

# R 4.6.1:  set.seed(100); rnorm(6)  and  set.seed(100); runif(4)
R_RNORM_100 = np.array([-0.502192350531456944, 0.131531165327302524, -0.078917089818869982,
                        0.886784809417844544, 0.116971270510841260, 0.318630087617032098])
R_RUNIF_100 = np.array([0.307766110869124532, 0.257672501029446721, 0.552322433330118656,
                        0.056383150396868587])


@contextlib.contextmanager
def _patched(mod, name, fn):
    old = getattr(mod, name)
    setattr(mod, name, fn)
    try:
        yield
    finally:
        setattr(mod, name, old)


def _fail(*args, **kwargs):
    raise np.linalg.LinAlgError("forced failure (test)")


def _fail_after(fn, ncalls):
    """``fn`` for the first ``ncalls`` calls, then a failure."""
    count = [0]

    def wrapped(*args, **kwargs):
        count[0] += 1
        if count[0] > ncalls:
            _fail()
        return fn(*args, **kwargs)
    return wrapped


def _expect_runtime_error(f, text):
    try:
        f()
    except RuntimeError as e:
        assert text in str(e), str(e)
    else:
        raise AssertionError("no RuntimeError")


# ------------------------------------------------------------ initial values ---
def test_r_random_numbers():
    assert np.array_equal(_core.r_runif(4, 100), R_RUNIF_100)
    assert np.allclose(_core.r_rnorm(6, 100), R_RNORM_100, rtol=0, atol=1e-15)
    assert _core.r_rnorm(0, 100).shape == (0,)


def test_default_beta_starts_are_rs_and_leave_numpy_rng_alone():
    ref = np.r_[0.0, 1.0 + R_RNORM_100[:5]]        # R: set.seed(100); c(0, rnorm(5, 1, 1))
    assert np.allclose(capcov.cap_beta0_v5(6), ref, rtol=0, atol=1e-15)
    assert np.allclose(capcov.coc_beta0_v3(6), ref, rtol=0, atol=1e-15)
    np.random.seed(7)
    a = np.random.rand(3)
    np.random.seed(7)
    capcov.cap_beta0_v5(6)
    capcov.coc_beta0_v3(4)
    b = np.random.rand(3)
    assert np.array_equal(a, b)
    for f in (hdcap._beta0_mat, coc._beta0_mat):
        B0 = f(6)
        assert B0.shape == (6, 2) and np.all(B0[:, 0] == 0) and np.allclose(B0[:, 1], ref)


def test_beta0_argument_checks():
    for f in (hdcap._beta0_mat, coc._beta0_mat):
        for bad in (np.zeros(2), np.array([0.0, np.nan, 1.0]), np.zeros((2, 2))):
            try:
                f(3, bad)
            except ValueError:
                pass
            else:
                raise AssertionError("bad beta0 accepted")
        assert f(3, [0, 1, 2]).shape == (3, 1)
        assert f(3, np.zeros((3, 4))).shape == (3, 4)


# -------------------------------------------------------------------- HDCAP ---
def test_hdcap_objective_is_the_function_the_fit_minimizes():
    Y, X, g, _ = _sim_cap()
    beta = np.array([0.1, 0.5, -0.2])
    Xb = X @ beta
    # without shrinkage: the negative log-likelihood, sample (MLE) covariances
    S, T = _core.sigma_cube(Y)
    ref = np.sum(T * (Xb + _core.scores(S, g) * np.exp(-Xb)))
    assert np.isclose(capcov.obj_func(Y, X, g, beta), ref, rtol=1e-12)
    # with shrinkage: per-subject cov_ls and the shrunk score of the beta fit
    Ssk, _ = _core.sigma_cube(Y, cov_shrinkage=True)
    sc = _core.scores(Ssk, g)
    w = hdcap._shrink_weights(sc, Xb, T, g @ g)
    ref_sk = np.sum(T * (Xb + (w["rho1"] * (g @ g) + w["rho2"] * sc) * np.exp(-Xb)))
    assert np.isclose(capcov.obj_func(Y, X, g, beta, cov_shrinkage=True), ref_sk, rtol=1e-12)


def test_cap_beta_minimizes_the_objective():
    Y, X, g, _ = _sim_cap()
    b = capcov.cap_beta(Y, X, g, cov_shrinkage=False)["beta"]
    f0 = capcov.obj_func(Y, X, g, b)
    for k in range(len(b)):
        e = np.zeros(len(b))
        e[k] = 1e-3
        assert capcov.obj_func(Y, X, g, b + e) > f0
        assert capcov.obj_func(Y, X, g, b - e) > f0


def test_cap_beta_several_starts_keep_the_converged_best():
    Y, X, g, _ = _sim_cap()
    good = capcov.cap_beta(Y, X, g, cov_shrinkage=False, beta0=np.zeros(3))
    bad = np.array([0.0, 400.0, -400.0])               # overflows: this start fails
    both = capcov.cap_beta(Y, X, g, cov_shrinkage=False,
                           beta0=np.column_stack([bad, np.zeros(3)]))
    assert good["convergence"] and both["convergence"]
    assert np.allclose(both["beta"], good["beta"], rtol=0, atol=1e-6)
    dflt = capcov.cap_beta(Y, X, g, cov_shrinkage=False)          # 0 and cap_beta0_v5
    assert np.allclose(dflt["beta"], good["beta"], rtol=0, atol=1e-6)
    _expect_runtime_error(lambda: capcov.cap_beta(Y, X, g, cov_shrinkage=False, beta0=bad),
                          "could not be estimated")


def test_cap_reg_beta0_default_is_zero_and_v5_start():
    Y, X, g_true, _ = _sim_cap()
    G0 = np.linalg.qr(np.random.RandomState(3).normal(size=(8, 4)))[0]
    r1 = capcov.cap_reg(Y, X, nD=2, gamma0_mat=G0)
    B0 = np.column_stack([np.zeros(3), capcov.cap_beta0_v5(3)])
    r2 = capcov.cap_reg(Y, X, nD=2, gamma0_mat=G0, beta0=B0)
    assert np.array_equal(r1["gamma"], r2["gamma"]) and np.array_equal(r1["beta"], r2["beta"])
    r0 = capcov.cap_reg(Y, X, nD=1, gamma0_mat=G0, beta0=np.zeros(3))
    assert abs(r0["gamma"][:, 0] @ g_true) > 0.99


def test_cap_reg_failures():
    Y, X, _, _ = _sim_cap()
    with _patched(hdcap, "_cap_d1", _fail):
        _expect_runtime_error(lambda: capcov.cap_reg(Y, X, nD=2), "cap_D1_opt")
    with _patched(hdcap, "_cap_dk", _fail):
        assert capcov.cap_reg(Y, X, nD=3)["nD"] == 1
        assert capcov.cap_reg(Y, X, stop_crt="DfD")["nD"] == 1


# ------------------------------------------------------------------ CAP-CoC ---
def test_coc_objective_uses_the_covariances_of_the_fit():
    d = _sim_coc()
    Y, X, W, g, t = d["Y"], d["X"], d["W"], d["gamma"], d["theta"]
    fit = capcov.coc_coef(Y, X, W, g, t)               # outcome shrinkage (default)
    Sx = coc._Sx_cube(X, True)
    Sy = coc._Sy_cube(Y, True)                         # per-subject cov_ls (start of the fit)
    ob = coc._obj_func(Sy, Sx, W, g, t, fit["alpha"], fit["beta"], Y=Y, cov_shrinkage_y=True)
    ref = coc._obj_scores(fit["score_y"], fit["score_x"], W, fit["alpha"], fit["beta"])
    assert np.isclose(ob, ref, rtol=1e-12)
    ob_v3 = coc._obj_func(Sy, Sx, W, g, t, fit["alpha"], fit["beta"])   # V3: cov_ls
    assert not np.isclose(ob_v3, ob, rtol=1e-6)


def test_coc_coef_several_starts_keep_the_converged_best():
    d = _sim_coc()
    Y, X, W, g, t = d["Y"], d["X"], d["W"], d["gamma"], d["theta"]
    singles = [capcov.coc_coef(Y, X, W, g, t, beta0=b) for b in coc._beta0_mat(3).T]
    obs = [coc._obj_scores(f["score_y"], f["score_x"], W, f["alpha"], f["beta"]) for f in singles]
    best = capcov.coc_coef(Y, X, W, g, t)              # both default starts
    k = int(np.argmin(obs))
    assert best["convergence"] and np.allclose(best["beta"], singles[k]["beta"])
    bad = np.array([0.0, 1e6, -1e6])
    both = capcov.coc_coef(Y, X, W, g, t, beta0=np.column_stack([bad, np.zeros(3)]))
    assert np.allclose(both["beta"], singles[0]["beta"], rtol=0, atol=1e-3)


def test_coc_reg_failures():
    d = _sim_coc()
    args = (d["Y"], d["X"], d["W"])
    kw = dict(burn_in=50, max_itr=300, ninitial=2)
    with _patched(coc, "_coc_d1", _fail):
        _expect_runtime_error(lambda: capcov.coc_reg(*args, nD=2, **kw), "COCReg.D1.opt")
    with _patched(coc, "_coc_dk", _fail):
        assert capcov.coc_reg(*args, nD=3, **kw)["gamma"].shape[1] == 1


# ------------------------------------------------- LCAP, mediation, clustering ---
def test_lcap_failures():
    d = _sim_lcap()
    with _patched(lcap, "_cap_d1", _fail):
        _expect_runtime_error(lambda: lcap.cap_reg(d["Y"], d["X"], nD=2, ninitial=2),
                              "cap_D1_opt")
    with _patched(lcap, "_cap_d1_opt", _fail_after(lcap._cap_d1_opt, 1)):
        fit = lcap.cap_reg(d["Y"], d["X"], nD=3, cov_shrinkage=False, ninitial=2)
    assert fit["gamma"].shape[1] == 1


def test_mediation_failures():
    d = _sim_med()
    with _patched(mediation, "cap_med_d1", _fail):
        _expect_runtime_error(lambda: mediation.cap_med(d["X"], d["M"], d["Y"], nD=2, ninitial=2),
                              "CAPMediation_D1_opt")
    with _patched(mediation, "_cap_med_d1_opt", _fail_after(mediation._cap_med_d1_opt, 1)):
        fit = mediation.cap_med(d["X"], d["M"], d["Y"], nD=3, ninitial=2)
    assert fit["theta"].shape[1] == 1


def test_clustering_failures():
    d = _sim_pcl(n=40, Ti=60)
    with _patched(clustering, "_cap_pcl_d1", _fail):
        _expect_runtime_error(lambda: clustering.cap_pcl(d["Y"], d["X"], d["W"], nD=2, ninitial=2),
                              "capPCL_D1_opt")
    with _patched(clustering, "_cap_pcl_dk", _fail):
        fit = clustering.cap_pcl(d["Y"], d["X"], d["W"], nD=3, ninitial=2)
    assert fit["gamma"].shape[1] == 1
