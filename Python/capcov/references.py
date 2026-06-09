"""Publication references for the CAP family of methods.

Source of truth: ``setting.md`` in the project root. Each entry gives the
citation, DOI, and URL for the method's paper. Methods without a publication yet
(MCAP, CAP-HDcov, CAP-clustering) are listed with ``None``.
"""
from __future__ import annotations

__all__ = ["REFERENCES", "references"]

REFERENCES = {
    "hdcap": [
        {
            "citation": ("Zhao, Y., Caffo, B., Luo, X., & Alzheimer's Disease "
                         "Neuroimaging Initiative (2021). Principal regression for "
                         "high dimensional covariance matrices. Electronic Journal "
                         "of Statistics, 15(2), 4192."),
            "doi": "10.1214/21-EJS1887",
            "url": "https://doi.org/10.1214/21-EJS1887",
        },
        {
            "citation": ("Zhao, Y., Wang, B., Mostofsky, S. H., Caffo, B. S., & "
                         "Luo, X. (2021). Covariate assisted principal regression "
                         "for covariance matrix outcomes. Biostatistics, 22(3), "
                         "629-645. [the classical CAP that HDCAP subsumes]"),
            "doi": "10.1093/biostatistics/kxz057",
            "url": "https://doi.org/10.1093/biostatistics/kxz057",
        },
    ],
    # CAP is obtained from HDCAP with cov_shrinkage=False -> same references.
    "cap": "hdcap",
    "coc": [
        {
            "citation": ("Zhao, Y., & Zhao, Y. (2025). Covariance-on-covariance "
                         "regression. Biometrics, 81(3), ujaf097."),
            "doi": "10.1093/biomtc/ujaf097",
            "url": "https://doi.org/10.1093/biomtc/ujaf097",
        },
    ],
    "lcap": [
        {
            "citation": ("Zhao, Y., Caffo, B. S., & Luo, X. (2024). Longitudinal "
                         "regression of covariance matrix outcomes. Biostatistics, "
                         "25(2), 385-401."),
            "doi": "10.1093/biostatistics/kxac045",
            "url": "https://doi.org/10.1093/biostatistics/kxac045",
        },
    ],
    "mediation": [
        {
            "citation": ("Xu, Y., & Zhao, Y. (2025). Mediation analysis with graph "
                         "mediator. Biostatistics, 26(1), kxaf004."),
            "doi": "10.1093/biostatistics/kxaf004",
            "url": "https://doi.org/10.1093/biostatistics/kxaf004",
        },
    ],
    "mcap": None,        # no publication yet
    "clustering": None,  # no publication yet
    "hdcov": None,       # no publication yet
}


def references(method=None):
    """Return (and pretty-print) the reference(s) for a method, or all methods.

    Parameters
    ----------
    method : str or None
        One of the keys of ``REFERENCES`` (e.g. ``"hdcap"``, ``"coc"``). If
        ``None``, prints every method's references and returns the full dict.
    """
    def _resolve(key):
        v = REFERENCES.get(key)
        return REFERENCES.get(v) if isinstance(v, str) else v

    if method is None:
        for key in REFERENCES:
            if isinstance(REFERENCES[key], str):
                continue
            refs = REFERENCES[key]
            print(f"[{key}]")
            if refs is None:
                print("  (no publication yet)")
            else:
                for r in refs:
                    print(f"  {r['citation']}")
                    print(f"  {r['url']}")
            print()
        return REFERENCES

    refs = _resolve(method)
    if refs is None:
        print(f"[{method}] (no publication yet)")
    else:
        for r in refs:
            print(r["citation"])
            print(r["url"])
    return refs
