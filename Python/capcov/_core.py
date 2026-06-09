"""Shared covariance-regression machinery for the CAP family.

These mirror the R helpers used across the CAP methods (cov.ls linear shrinkage,
per-subject covariances, projected scores, weighted accumulation, the
generalized smallest-eigenvector solve, and the deviation-from-diagonality).
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

__all__ = [
    "as_matrix_list", "cov_ls", "sigma_cube", "scores", "accum",
    "gamma_solve", "diag_level",
]


def as_matrix_list(Y):
    """Coerce a list/sequence of array-likes to a list of 2-D float arrays."""
    return [np.asarray(Yi, dtype=float).reshape(np.asarray(Yi).shape[0], -1) for Yi in Y]


def cov_ls(X):
    """Ledoit--Wolf style linear shrinkage of one ``n x p`` data matrix.

    Matches the R ``cov.ls``: shrink the (MLE) sample covariance toward a scaled
    identity, ``b2*m/d2 * I + a2/d2 * S``.
    """
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    Xc = X - X.mean(axis=0, keepdims=True)
    S = (Xc.T @ Xc) / n                      # cov * (n-1)/n
    m = np.trace(S) / p                      # norm.F.std(S, I)
    A = S - m * np.eye(p)
    d2 = np.sum(A * A) / p                    # norm.F.std(S - m I)^2
    normS2 = np.sum(S * S)
    xx = np.einsum("ij,ij->i", Xc, Xc)       # x_i' x_i
    xSx = np.einsum("ij,jk,ik->i", Xc, S, Xc)  # x_i' S x_i
    b2_bar = np.mean((xx * xx - 2.0 * xSx + normS2) / p) / n
    b2 = min(b2_bar, d2)
    a2 = d2 - b2
    return b2 * m * np.eye(p) / d2 + a2 * S / d2


def sigma_cube(Y, cov_shrinkage=False, second_moment=False):
    """Per-subject covariance array of shape ``(p, p, n)`` plus sample sizes.

    - ``cov_shrinkage`` : use ``cov_ls`` linear shrinkage instead of the sample
      covariance (the HDCAP option).
    - ``second_moment`` : use ``Y_i' Y_i / T_i`` (no centering) instead of the
      centered covariance (used by CAP-clustering).
    """
    Y = as_matrix_list(Y)
    n = len(Y)
    p = Y[0].shape[1]
    cube = np.empty((p, p, n))
    Tvec = np.empty(n)
    for i, Yi in enumerate(Y):
        Ti = Yi.shape[0]
        Tvec[i] = Ti
        if cov_shrinkage:
            cube[:, :, i] = cov_ls(Yi)
        elif second_moment:
            cube[:, :, i] = (Yi.T @ Yi) / Ti
        else:
            Yc = Yi - Yi.mean(axis=0, keepdims=True)
            cube[:, :, i] = (Yc.T @ Yc) / Ti
    return cube, Tvec


def scores(Sigma, v):
    """Projected variances ``v' Sigma_i v`` for all ``i`` (``Sigma`` is p x p x n)."""
    v = np.asarray(v, dtype=float).ravel()
    return np.einsum("i,ijk,j->k", v, Sigma, v)


def accum(Sigma, w):
    """Weighted accumulation ``sum_i w_i Sigma_i`` (``Sigma`` is p x p x n)."""
    w = np.asarray(w, dtype=float).ravel()
    return Sigma @ w


def gamma_solve(A, H):
    """Smallest generalized eigenvector: minimize ``g' A g`` s.t. ``g' H g = 1``.

    Mirrors the R ``gamma.solve`` / ``eigen.solve``. ``scipy.linalg.eigh(A, H)``
    returns ascending generalized eigenpairs with ``V' H V = I``; the first
    column is the minimizer.
    """
    A = np.asarray(A, dtype=float)
    H = np.asarray(H, dtype=float)
    A = 0.5 * (A + A.T)
    H = 0.5 * (H + H.T)
    vals, vecs = eigh(A, H)
    return vecs[:, 0]


def diag_level(Y, Phi):
    """Deviation-from-diagonality (DfD) for the first ``k`` directions.

    Returns ``(avg_level, sub_level)`` matching the R ``diag.level``:
    ``avg.level[k]`` is the size-weighted geometric mean across subjects of
    ``det(diag(M)) / det(M)`` where ``M = Phi[:, :k]' cov(Y_i) Phi[:, :k]``.
    """
    Y = as_matrix_list(Y)
    Phi = np.asarray(Phi, dtype=float)
    if Phi.ndim == 1:
        Phi = Phi[:, None]
    n = len(Y)
    ps = Phi.shape[1]
    Tvec = np.array([Yi.shape[0] for Yi in Y], dtype=float)
    dl = np.ones((n, ps))
    for i, Yi in enumerate(Y):
        Yc = Yi - Yi.mean(axis=0, keepdims=True)
        cov_i = (Yc.T @ Yc) / (Yi.shape[0] - 1)   # R cov() uses (T-1)
        for j in range(1, ps):
            P = Phi[:, : j + 1]
            M = P.T @ cov_i @ P
            dl[i, j] = np.prod(np.diag(M)) / np.linalg.det(M)
    avg = np.array([
        np.prod(dl[:, j] ** (Tvec / Tvec.sum())) for j in range(ps)
    ])
    return avg, dl
