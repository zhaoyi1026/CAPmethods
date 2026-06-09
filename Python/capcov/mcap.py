"""MCAP -- Multilevel CAP (cluster-varying projection, vMF)

Planned port of `LCAP_gamma-var/V4`. Not yet implemented in Python; use the R
implementation (`LCAP_gamma-var/V4`, RcppArmadillo-accelerated) meanwhile. See the package
README for the porting roadmap.
"""
from __future__ import annotations


def mcap_reg(*args, **kwargs):
    raise NotImplementedError(
        "MCAP -- Multilevel CAP (cluster-varying projection, vMF) is not yet ported to Python. "
        "Use the R implementation LCAP_gamma-var/V4 for now."
    )
