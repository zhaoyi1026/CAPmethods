"""Each built-in example generator returns the right shape and recovers truth.

Kept in sync with the R ``CAPmethods`` ``*_example()`` generators (manuscript /
Shiny-demo settings). Tests use light fit settings (small ``ninitial``, single
direction where the second is slow) — fuller recovery is shown in the package
README / the R vignette.
"""
import numpy as np
import capcov
from capcov import examples as ex


def test_references_cover_implemented_methods():
    REF = capcov.REFERENCES
    for m in ["hdcap", "coc", "lcap", "mediation"]:
        refs = REF[m]
        assert isinstance(refs, list) and refs and "doi" in refs[0]
    assert REF["cap"] == "hdcap"          # CAP shares HDCAP's refs
    assert REF["clustering"] is None       # no publication yet


def test_hdcap_example_recovers_two_directions():
    d = ex.hdcap_example()
    assert d["X"].shape == (100, 2) and len(d["Y"]) == 100 and d["Y"][0].shape[1] == 20
    assert d["truth"]["gamma"].shape == (20, 2) and d["truth"]["beta"].shape == (2, 2)
    f = capcov.cap_reg(d["Y"], d["X"], stop_crt="nD", nD=2, cov_shrinkage=True, ninitial=5)
    G = d["truth"]["gamma"]
    for k in range(2):                                  # each found dir is in the true 2-D span
        cos = max(abs(f["gamma"][:, k] @ G[:, 0]), abs(f["gamma"][:, k] @ G[:, 1]))
        assert cos > 0.9
    assert min(abs(f["beta"][1, :])) > 0.6             # group effect (|±1|) recovered


def test_lcap_example_recovers_leading_direction():
    d = ex.lcap_example()
    assert len(d["Y"]) == 100 and d["Y"][0][0].shape[1] == 20
    assert d["truth"]["gamma"].shape == (20, 2) and d["truth"]["beta"].shape == (3, 2)
    f = capcov.lcap.cap_reg(d["Y"], d["X"], stop_crt="nD", nD=1,
                            cov_shrinkage=True, ninitial=3)
    G = d["truth"]["gamma"]
    assert max(abs(f["gamma"][:, 0] @ G[:, 0]), abs(f["gamma"][:, 0] @ G[:, 1])) > 0.9


def test_coc_example_recovers():
    d = ex.coc_example()
    assert len(d["Y"]) == 150 and d["X"][0].shape[1] == 10 and d["Y"][0].shape[1] == 5
    rs = d["truth"]["recovering_settings"]
    f = capcov.coc_reg(d["Y"], d["X"], d["W"], stop_crt="nD", nD=1,
                       Hy=rs["Hy"], Hx=rs["Hx"], burn_in=rs["burn_in"], ninitial=5)
    g = f["gamma"][:, 0]; t = f["theta"][:, 0]
    gcos = max(abs(g @ d["truth"]["gamma"][:, 0]), abs(g @ d["truth"]["gamma"][:, 1]))
    tcos = max(abs(t @ d["truth"]["theta"][:, 0]), abs(t @ d["truth"]["theta"][:, 1]))
    assert gcos > 0.9 and tcos > 0.9


def test_mediation_example_recovers():
    d = ex.mediation_example()
    assert d["X"].shape == (100, 1) and len(d["M"]) == 100 and d["M"][0].shape == (150, 10)
    assert abs(d["truth"]["IE"] - 0.56) < 1e-9
    f = capcov.cap_med(d["X"], d["M"], d["Y"], stop_crt="nD", nD=1, ninitial=5)
    assert abs(f["theta"][:, 0] @ d["truth"]["theta"]) > 0.95   # one true direction
    assert abs(float(f["alpha"][0, 0]) - 0.8) < 0.2             # a-path recovered


def test_clustering_example_recovers_separated_component():
    d = ex.clustering_example()
    assert len(d["Y"]) == 100 and d["Y"][0].shape[1] == 50
    assert d["truth"]["gamma"].shape == (50, 2) and d["truth"]["cluster"].shape == (100, 2)
    g_D2 = d["truth"]["gamma"][:, 0]; cl_D2 = d["truth"]["cluster"][:, 0] - 1
    cf = capcov.clustering.cap_pcl_coef(d["Y"], d["X"], d["W"], g_D2, ncluster=2, nstart=8)
    agree = max(np.mean((cf["class"] - 1) == cl_D2), np.mean((cf["class"] - 1) == (1 - cl_D2)))
    assert agree > 0.9
