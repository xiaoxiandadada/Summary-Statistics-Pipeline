#!/usr/bin/env python3
"""
Lightweight Python accelerator for LAVA‑Knock pipeline.

Functions are designed to be imported from R via:
  reticulate::source_python('accelerator.py')

Provided capabilities:
- fast_correlation_matrix: Numba/NumPy accelerated correlation
- ld_pruning: greedy LD pruning on absolute correlation
- preprocess_genotype: QC and filter genotype matrix
- fast_svd_decomposition: truncated SVD by variance threshold
- compute_window_statistics: per-window univariate/bivariate stats
- knockoff_filter_fast: fast knockoff filter (W, threshold, selection)
- parallel_window_analysis: parallel per-window processing (optional)
"""

from typing import Optional, List, Dict, Tuple

try:
    import numpy as np
except Exception as e:  # pragma: no cover
    np = None

# Try optional numba JIT; fallback to numpy implementation if unavailable.
try:  # pragma: no cover
    from numba import jit, prange
    _HAVE_NUMBA = True
except Exception:
    _HAVE_NUMBA = False
    def jit(*args, **kwargs):  # type: ignore
        def _wrap(fn):
            return fn
        return _wrap
    def prange(*args, **kwargs):  # type: ignore
        return range(*args)


def _corr_numpy(X: "np.ndarray") -> "np.ndarray":
    # rows=samples, cols=features; ensure float64
    X = np.asarray(X, dtype=np.float64)
    # Use column-wise corrcoef
    return np.corrcoef(X, rowvar=False)

@jit(nopython=True, parallel=True)
def _corr_numba(X):  # type: ignore
    n, p = X.shape
    corr = np.zeros((p, p))
    # standardize columns
    for j in prange(p):
        mean_j = np.mean(X[:, j])
        std_j = np.std(X[:, j])
        if std_j > 0.0:
            X[:, j] = (X[:, j] - mean_j) / std_j
        else:
            X[:, j] = 0.0
    for i in prange(p):
        for j in range(i, p):
            v = 0.0
            for k in range(n):
                v += X[k, i] * X[k, j]
            v = v / (n - 1)
            corr[i, j] = v
            corr[j, i] = v
    return corr


def fast_correlation_matrix(X):
    """Compute correlation matrix for a samples×SNPs matrix.

    If Numba is available, uses a parallel JIT; otherwise falls back to numpy.
    Returns a 2D numpy array (float64).
    """
    if np is None:
        raise RuntimeError("numpy is required for accelerator")
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError("X must be a 2D array")
    if _HAVE_NUMBA:
        # make a copy because numba version modifies X in-place during standardization
        return _corr_numba(X.copy())
    return _corr_numpy(X)


def ld_pruning(corr_matrix, threshold: float = 0.75):
    """Greedy LD pruning on absolute correlation matrix.

    Keeps a subset of indices such that for each kept i, |corr[i, kept]| <= threshold
    in a single-pass greedy order. Returns a numpy array of 0-based kept indices.
    """
    if np is None:
        raise RuntimeError("numpy is required for accelerator")
    C = np.asarray(corr_matrix, dtype=np.float64)
    p = C.shape[0]
    kept = []
    for i in range(p):
        keep = True
        for j in kept:
            if abs(C[i, j]) > threshold:
                keep = False
                break
        if keep:
            kept.append(i)
    return np.array(kept, dtype=np.int64)


def available() -> bool:
    """Quick probe for availability from R."""
    return np is not None


# ----------------------------
# Additional utilities
# ----------------------------

def preprocess_genotype(G: "np.ndarray") -> "np.ndarray":
    """QC filter genotype matrix (samples×SNPs).

    - Clamp to [0,2]; set invalid to 0
    - Filter SNPs with maf>0, mac>=25, var>0
    Returns filtered copy of G (not inplace).
    """
    if np is None:
        raise RuntimeError("numpy is required for accelerator")
    X = np.asarray(G, dtype=np.float64)
    X = np.where((X < 0) | (X > 2) | ~np.isfinite(X), 0.0, X)
    maf = X.mean(axis=0) / 2.0
    mac = X.sum(axis=0)
    var = X.var(axis=0)
    mask = (maf > 0.0) & (mac >= 25.0) & (var > 0.0) & np.isfinite(maf)
    if mask.sum() < 2:
        return X[:, :0]
    return X[:, mask]


def fast_svd_decomposition(X: "np.ndarray", variance_threshold: float = 0.99) -> Tuple["np.ndarray","np.ndarray","np.ndarray"]:
    """Compute truncated SVD such that cumulated variance >= threshold.

    Returns (U_k, s_k, Vt_k).
    """
    if np is None:
        raise RuntimeError("numpy is required for accelerator")
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    if s.size == 0:
        return U, s, Vt
    ve = (s ** 2) / float((s ** 2).sum())
    cum = np.cumsum(ve)
    k = int(np.argmax(cum >= variance_threshold)) + 1
    return U[:, :k], s[:k], Vt[:k, :]


def compute_bivariate_pvalue(omega: "np.ndarray", sigma2: "np.ndarray", K: int, n_iter: int = 1000) -> float:
    """Approximate p-value for local rg using normal approximation.

    NOTE: Simplified; for production use replace with accurate Wishart-based test.
    """
    if omega.shape[0] < 2 or np.any(np.diag(omega) <= 0):
        return float("nan")
    obs_rg = omega[0, 1] / np.sqrt(omega[0, 0] * omega[1, 1])
    try:
        se_rg = np.sqrt(sigma2[0] * sigma2[1]) / (K * np.sqrt(omega[0, 0] * omega[1, 1]))
        z = obs_rg / se_rg if se_rg > 0 else 0.0
        # two-sided normal
        from math import erf, sqrt
        # 1 - Phi(|z|) = 0.5 * erfc(|z|/sqrt(2)); use approximation via erf
        # Here convert to p = 2*(1-Phi(|z|))
        # Phi(z) ~ 0.5*(1+erf(z/sqrt(2)))
        p_one = 0.5 * (1.0 - erf(abs(z) / sqrt(2.0)))
        return 2.0 * p_one
    except Exception:
        return float("nan")


def compute_window_statistics(genotype_matrix: "np.ndarray", zscore_matrix: "np.ndarray", n_samples: int) -> Dict[str, float]:
    """Compute per-window stats: univariate h2/p and bivariate rg/p.

    Inputs:
      genotype_matrix: samples×SNPs
      zscore_matrix: SNPs×2 (pheno1, pheno2) or SNPs×k
    """
    if np is None:
        raise RuntimeError("numpy is required for accelerator")

    # QC
    G = preprocess_genotype(genotype_matrix)
    if G.shape[1] < 2:
        return {"h2_pheno1": np.nan, "h2_pheno2": np.nan, "p_uni_pheno1": np.nan, "p_uni_pheno2": np.nan, "rg": np.nan, "bivar_p": np.nan, "K": 0}

    # Correlation/SVD basis
    # Note: to be consistent with R pipeline, use genotype-based PCA basis
    X = (G - G.mean(axis=0, keepdims=True)) / (G.std(axis=0, keepdims=True) + 1e-8)
    U, s, Vt = fast_svd_decomposition(X.T / np.sqrt(max(1, n_samples - 1)), variance_threshold=0.99)
    K = len(s)
    if K == 0:
        return {"h2_pheno1": np.nan, "h2_pheno2": np.nan, "p_uni_pheno1": np.nan, "p_uni_pheno2": np.nan, "rg": np.nan, "bivar_p": np.nan, "K": 0}

    Z = np.asarray(zscore_matrix, dtype=np.float64)
    if Z.ndim == 1:
        Z = Z.reshape((-1, 1))
    # Simple projection-based approximation of LAVA quantities
    Q = Vt.T
    lam = s ** 2
    try:
        alpha = Q @ np.diag(1.0 / lam) @ Q.T @ Z
    except Exception:
        # numerical fallback
        alpha = Q @ np.linalg.pinv(np.diag(lam)) @ Q.T @ Z
    delta = (np.sqrt(lam)[:, None] * (Q.T @ alpha))
    R2 = np.diag(Z.T @ alpha)
    eta2 = (n_samples - 1) / max(1, (n_samples - K - 1)) * (1 - R2)
    sigma2 = eta2 / max(1, (n_samples - 1))
    h2 = 1 - eta2
    T_uni = np.diag(delta.T @ delta) / (sigma2 * max(1, K))
    try:
        from scipy.stats import f as fdist
        p_uni = 1 - fdist.cdf(T_uni, K, max(1, n_samples - K - 1))
    except Exception:
        p_uni = np.full_like(T_uni, np.nan, dtype=float)

    if Z.shape[1] >= 2:
        omega = (delta.T @ delta) / max(1, K) - np.diag(sigma2)
        try:
            rg = omega[0, 1] / np.sqrt(omega[0, 0] * omega[1, 1])
        except Exception:
            rg = float("nan")
        bivar_p = compute_bivariate_pvalue(omega, sigma2, K)
    else:
        rg, bivar_p = float("nan"), float("nan")

    return {
        "h2_pheno1": float(h2[0]) if h2.size > 0 else float("nan"),
        "h2_pheno2": float(h2[1]) if h2.size > 1 else float("nan"),
        "p_uni_pheno1": float(p_uni[0]) if p_uni.size > 0 else float("nan"),
        "p_uni_pheno2": float(p_uni[1]) if p_uni.size > 1 else float("nan"),
        "rg": float(rg),
        "bivar_p": float(bivar_p) if bivar_p == bivar_p else float("nan"),
        "K": int(K),
    }


def knockoff_filter_fast(p_original: "np.ndarray", p_knockoffs: "np.ndarray", fdr: float = 0.1, M: int = 5) -> Dict[str, "np.ndarray"]:
    """Fast knockoff filter aggregation.

    Inputs:
      p_original: shape (n,)
      p_knockoffs: shape (n, M)
    Returns dict with W, threshold, significant (bool), n_significant.
    """
    if np is None:
        raise RuntimeError("numpy is required for accelerator")
    p0 = np.asarray(p_original, dtype=np.float64).reshape((-1,))
    pko = np.asarray(p_knockoffs, dtype=np.float64)
    if pko.ndim == 1:
        pko = pko.reshape((-1, 1))
    T0 = -np.log10(np.clip(p0, 1e-300, 1.0))
    T_ko = -np.log10(np.clip(pko, 1e-300, 1.0))
    T = np.column_stack([T0, T_ko])
    T[np.isnan(T)] = 0.0
    if M > 1:
        T_ko_median = np.median(T_ko, axis=1)
        T_ko_max = np.max(T_ko, axis=1)
        W = (T0 - T_ko_median) * (T0 >= T_ko_max)
    else:
        W = T0 - T_ko[:, 0]
    kappa = np.argmax(T, axis=1)
    tau = np.max(T, axis=1) - np.median(np.delete(T, np.argmax(T, axis=1), axis=1), axis=1)
    thr = compute_knockoff_threshold(kappa, tau, M, fdr)
    sig = W >= thr
    return {"W": W, "threshold": float(thr), "significant": sig.astype(bool), "n_significant": int(sig.sum())}


def compute_knockoff_threshold(kappa: "np.ndarray", tau: "np.ndarray", M: int, fdr: float, rej_bound: int = 20000) -> float:
    order = np.argsort(-tau)
    c0 = (kappa[order] == 0)
    ratios = []
    temp0 = 0
    L = min(len(order), rej_bound)
    for i in range(L):
        temp0 += int(c0[i])
        temp1 = (i + 1) - temp0
        ratios.append((1.0 / M + (1.0 / M) * temp1) / max(1, temp0))
    arr = np.asarray(ratios, dtype=np.float64)
    idx = np.where(arr <= fdr)[0]
    if idx.size > 0:
        return float(tau[order[idx[-1]]])
    return float("inf")


def parallel_window_analysis(genotype_data: "np.ndarray", zscore_df, windows_df, window_size: int = 100000, n_samples: int = 20000) -> List[Dict]:
    """Parallel window analysis. Windows dataframe must have columns ['start','end'].
    zscore_df must have at least columns ['pos','zscore_pheno1','zscore_pheno2'] and be index-aligned to genotype columns.
    """
    try:
        import pandas as pd  # type: ignore
    except Exception:  # pragma: no cover
        raise RuntimeError("pandas is required for parallel_window_analysis")
    from concurrent.futures import ProcessPoolExecutor

    def _analyze(args):
        widx, wstart, wend = args
        mask = (zscore_df['pos'] >= wstart) & (zscore_df['pos'] <= wend)
        idx = np.where(mask.values)[0]
        if idx.size < 5:
            return None
        Gwin = genotype_data[:, idx]
        Zwin = zscore_df.loc[mask, ['zscore_pheno1', 'zscore_pheno2']].values
        res = compute_window_statistics(Gwin, Zwin, n_samples)
        res.update({"window_idx": int(widx), "window_start": int(wstart), "window_end": int(wend), "n_variants": int(idx.size)})
        return res

    args_iter = [(i, int(w['start']), int(w['end'])) for i, w in windows_df.iterrows()]
    with ProcessPoolExecutor() as ex:
        out = list(ex.map(_analyze, args_iter))
    return [r for r in out if r is not None]


# Note: CLI removed to avoid interference when imported via reticulate (which sets __name__=='__main__').
