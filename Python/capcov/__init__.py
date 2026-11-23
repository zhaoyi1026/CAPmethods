"""capcov -- Covariate Assisted Principal (CAP) family of covariance-regression
methods, ported from the R reference implementations.

CAP itself is not a separate module: it is ``hdcap.cap_reg(..., cov_shrinkage=False)``.

The top-level ``cap_reg`` / ``cap_beta`` / ``cap_beta_boot`` are HDCAP's. LCAP and
CAP-clustering reuse those names for their own models, so they are exposed as
submodules: ``capcov.lcap.cap_reg(...)``, ``capcov.clustering.cap_pcl(...)``.

0.2.0 follows CAPmethods 0.2.0 (HDCAP V9, CAP-CoC V4): the initial values of a
direction are compared by the function the fit minimizes, beta has two initial
values by default (``beta0``), and a direction that cannot be estimated gives a
clear error (first direction) or ends the search (later directions).

Method status (see README):
    hdcap      -- HDCAP / CAP            : implemented, verified vs R
    coc        -- CAP-CoC               : implemented, verified vs R
    lcap       -- LCAP (gamma-invariant): implemented, verified vs R
    mcap       -- MCAP (gamma-varying)  : planned
    mediation  -- CAP-mediation         : implemented (approx.), verified vs R
    clustering -- CAP-clustering        : implemented (approx.), truth-verified
    hdcov      -- CAP-HDcov (HCAP)      : planned
"""
from __future__ import annotations

from . import _core
from . import lcap
from . import clustering
from . import mediation
from . import examples
from . import references as _references
from .hdcap import cap_reg, cap_beta, cap_beta_boot, obj_func, cap_beta0_v5
from .coc import coc_reg, coc_coef, coc_coef_asmp, coc_coef_boot, coc_beta0_v3
from .mediation import cap_med, cap_med_coef, cap_med_d1
from .references import REFERENCES, references

__all__ = [
    "cap_reg", "cap_beta", "cap_beta_boot", "obj_func", "cap_beta0_v5",
    "coc_reg", "coc_coef", "coc_coef_asmp", "coc_coef_boot", "coc_beta0_v3",
    "cap_med", "cap_med_coef", "cap_med_d1",
    "lcap", "clustering", "mediation", "examples",
    "REFERENCES", "references", "_core",
]
__version__ = "0.2.0"
