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

from typing import Optional, List, Dict, Tuple, Sequence
import os
import csv
import warnings
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor

try:
    import numpy as np
except Exception as e:  # pragma: no cover
    np = None

try:  # pragma: no cover
    import pandas as pd
except Exception:
    pd = None

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
# Helpers for pandas-based preprocessing
# ----------------------------

def _require_pandas():
    if pd is None:
        raise RuntimeError("pandas is required for Python-side preprocessing")


def _canon_chr(series):
    _require_pandas()
    s = series.astype(str).str.upper().str.replace('^CHR', '', regex=True)
    s = s.replace({'X': '23', 'Y': '24', 'M': '25', 'MT': '25', 'MITO': '25'})
    return pd.to_numeric(s, errors='coerce')


def _canon_pos(series):
    _require_pandas()
    return pd.to_numeric(series, errors='coerce')


def _read_table_auto(path):
    _require_pandas()
    return pd.read_csv(path, sep=None, engine='python')


# ----------------------------
# LD block partitioning helpers
# ----------------------------

BLOCK_FILES = {
    "GRCH37": "LAVA_s2500_m25_f1_w200.blocks",
    "GRCH38": "deCODE_EUR_LD_blocks.bed"
}

GENE_FILES = {
    "GRCH37": "coding.genes.TSS.hg19.tsv",
    "GRCH38": "coding.genes.hg38.tsv"
}


def _resolve_block_path(coord_version: str, base_dir: str = ".") -> str:
    candidate = coord_version.strip()
    if os.path.isfile(candidate):
        return os.path.abspath(candidate)
    key = candidate.upper()
    if key in BLOCK_FILES:
        filename = BLOCK_FILES[key]
        candidate_paths = [filename, os.path.join(base_dir, filename)]
        path = next((p for p in candidate_paths if os.path.exists(p)), None)
        if path is None:
            raise FileNotFoundError(f"Cannot locate LD block file: {filename}")
        return os.path.abspath(path)
    raise FileNotFoundError(f"Unsupported coord version or path: {coord_version}")


@lru_cache(maxsize=16)
def _load_ld_blocks_cached(path: str) -> Dict[int, List[Tuple[float, float]]]:
    per_chr: Dict[int, List[Tuple[float, float]]] = {}
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chr_val = row.get("chr") or row.get("CHR") or row.get("chrom")
            if chr_val is None:
                continue
            chr_val = chr_val.strip()
            if chr_val.lower().startswith("chr"):
                chr_val = chr_val[3:]
            try:
                chr_int = int(chr_val)
            except ValueError:
                continue
            start = row.get("start") or row.get("Start")
            end = row.get("stop") or row.get("end") or row.get("Stop")
            if start is None or end is None:
                continue
            try:
                start_f = float(start)
                end_f = float(end)
            except ValueError:
                continue
            per_chr.setdefault(chr_int, []).append((min(start_f, end_f), max(start_f, end_f)))
    return per_chr


def partition_ld_blocks(
    chromosomes: List[int],
    positions: List[float],
    coord_version: str,
    fallback: bool = True,
    workers: Optional[int] = None,
    base_dir: str = "."
) -> Dict[str, List[Dict[str, object]]]:
    if np is None:
        raise RuntimeError("numpy is required for partition_ld_blocks")
    block_path = _resolve_block_path(coord_version, base_dir)
    blocks = _load_ld_blocks_cached(block_path)
    chrom_arr = np.asarray(chromosomes, dtype=np.int64)
    pos_arr = np.asarray(positions, dtype=np.float64)
    valid_mask = np.isfinite(chrom_arr) & np.isfinite(pos_arr)
    valid_indices = np.nonzero(valid_mask)[0]
    if valid_indices.size == 0:
        return {"chunks": [], "fallback": []}
    chrom_arr = chrom_arr[valid_indices]
    pos_arr = pos_arr[valid_indices]
    orig_idx = valid_indices
    unique_chr = np.unique(chrom_arr)
    max_workers = workers or min(len(unique_chr), os.cpu_count() or 1)

    def process_chr(chr_val: int):
        chr_mask = chrom_arr == chr_val
        chr_indices = orig_idx[chr_mask]
        chr_positions = pos_arr[chr_mask]
        if chr_indices.size == 0:
            return [], []
        assigned = np.zeros(chr_indices.shape[0], dtype=bool)
        chunk_list: List[Dict[str, object]] = []
        fallback_list: List[Dict[str, object]] = []
        chr_blocks = blocks.get(int(chr_val), [])
        for block_start, block_end in chr_blocks:
            block_mask = (chr_positions >= block_start) & (chr_positions <= block_end)
            if int(block_mask.sum()) < 5:
                continue
            block_indices = chr_indices[block_mask].tolist()
            chunk_list.append({
                "chr": int(chr_val),
                "start": float(block_start),
                "end": float(block_end),
                "indices": block_indices,
                "source": "ld_block"
            })
            assigned |= block_mask
        if fallback:
            leftover_mask = ~assigned
            if int(leftover_mask.sum()) >= 5:
                leftover_indices = chr_indices[leftover_mask].tolist()
                leftover_positions = chr_positions[leftover_mask]
                fallback_list.append({
                    "chr": int(chr_val),
                    "start": float(np.min(leftover_positions)),
                    "end": float(np.max(leftover_positions)),
                    "indices": leftover_indices,
                    "source": "fallback"
                })
        return chunk_list, fallback_list

    results_chunks: List[Dict[str, object]] = []
    results_fallback: List[Dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=max_workers or 1) as executor:
        for chunk_list, fallback_list in executor.map(process_chr, unique_chr.tolist()):
            if chunk_list:
                results_chunks.extend(chunk_list)
            if fallback_list:
                results_fallback.extend(fallback_list)
    return {"chunks": results_chunks, "fallback": results_fallback}


# ----------------------------
# Genotype loaders (R wraps these via reticulate)
# ----------------------------

_GENO_LUT = None


def _build_geno_lut():
    global _GENO_LUT
    if _GENO_LUT is not None:
        return _GENO_LUT
    lut = np.empty((256, 4), dtype=np.float64)
    lut.fill(np.nan)
    for byte_val in range(256):
        for shift in range(4):
            code = (byte_val >> (2 * shift)) & 0b11
            if code == 0:
                lut[byte_val, shift] = 0.0
            elif code == 1:
                lut[byte_val, shift] = np.nan  # missing
            elif code == 2:
                lut[byte_val, shift] = 1.0
            else:
                lut[byte_val, shift] = 2.0
    _GENO_LUT = lut
    return lut


def _read_plink_bed_subset(path: str, n_samples: int, selected_idx: "np.ndarray") -> "np.ndarray":
    lut = _build_geno_lut()
    selected_idx = np.asarray(selected_idx, dtype=np.int64)
    if selected_idx.ndim != 1:
        raise ValueError("selected_idx must be 1D array")
    if selected_idx.size == 0:
        return np.zeros((n_samples, 0), dtype=np.float64)
    order = np.argsort(selected_idx)
    selected_idx = selected_idx[order]
    bytes_per_snp = (n_samples + 3) // 4
    with open(path, "rb") as handle:
        header = handle.read(3)
        if len(header) != 3 or header[0] != 0x6C or header[1] != 0x1B:
            raise ValueError("Invalid PLINK .bed header (magic bytes mismatch)")
        if header[2] != 0x01:
            raise ValueError("PLINK .bed must be in SNP-major mode (header byte != 1)")
        result = np.empty((selected_idx.size, n_samples), dtype=np.float64)
        prev_idx = -1
        for out_pos, snp_idx in enumerate(selected_idx):
            if snp_idx <= prev_idx:
                raise ValueError("selected indices must be strictly increasing")
            skip = snp_idx - prev_idx - 1
            if skip > 0:
                handle.seek(skip * bytes_per_snp, os.SEEK_CUR)
            raw = np.fromfile(handle, dtype=np.uint8, count=bytes_per_snp)
            if raw.size != bytes_per_snp:
                raise ValueError("Unexpected EOF while reading PLINK bed")
            decoded = lut[raw].reshape(-1)[:n_samples]
            result[out_pos, :] = decoded
            prev_idx = snp_idx
    return result.T[:, order]


def load_genotype_rds(path: str):  # pragma: no cover (handled in R)
    raise NotImplementedError("RDS genotype loading is handled in R (requires R serialization)")


def load_genotype_csv(path: str, fmt: str = "samples_by_snps") -> "np.ndarray":
    if np is None:
        raise RuntimeError("numpy is required for load_genotype_csv")
    import pandas as pd
    df = pd.read_csv(path)
    mat = df.to_numpy(dtype=float)
    if fmt == "snps_by_samples":
        mat = mat.T
    return mat


def load_genotype_plink(path_prefix: str, snp_ids=None, chrpos_ids=None):
    if np is None or pd is None:
        raise RuntimeError("numpy and pandas are required for load_genotype_plink")

    bed_path = f"{path_prefix}.bed"
    bim_path = f"{path_prefix}.bim"
    fam_path = f"{path_prefix}.fam"
    for file_path in (bed_path, bim_path, fam_path):
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Missing PLINK file: {file_path}")

    bim_cols = ['CHR', 'RSID', 'CM', 'POS', 'A1', 'A2']
    bim_df = pd.read_csv(bim_path, sep=r'\s+', names=bim_cols)
    bim_df['CHR'] = _canon_chr(bim_df['CHR'])
    bim_df['POS'] = _canon_pos(bim_df['POS'])
    bim_df['CHRPOS'] = bim_df['CHR'].astype('Int64').astype(str) + ':' + bim_df['POS'].astype('Int64').astype(str)

    fam_df = pd.read_csv(fam_path, sep=r'\s+', header=None)
    n_samples = fam_df.shape[0]
    n_snps_total = bim_df.shape[0]

    select_mask = np.zeros(n_snps_total, dtype=bool)
    rsid_series = bim_df['RSID'].fillna('')
    chrpos_series = bim_df['CHRPOS']

    messages: List[str] = []

    if snp_ids:
        snp_list = [s for s in snp_ids if isinstance(s, str) and len(s) > 0]
        if snp_list:
            snp_set = set(snp_list)
            matched = rsid_series.isin(snp_set).to_numpy()
            select_mask |= matched
            missing_rsids = list(snp_set - set(rsid_series[matched]))
            if missing_rsids:
                messages.append("rsid missing: " + ", ".join(missing_rsids[:10]))

    if chrpos_ids:
        chrpos_list = [s for s in chrpos_ids if isinstance(s, str) and len(s) > 0]
        if chrpos_list:
            chrpos_set = set(chrpos_list)
            matched = chrpos_series.isin(chrpos_set).to_numpy()
            select_mask |= matched
            missing_chrpos = list(chrpos_set - set(chrpos_series[matched]))
            if missing_chrpos:
                messages.append("chr:pos missing: " + ", ".join(missing_chrpos[:10]))

    if not select_mask.any():
        select_mask = np.ones(n_snps_total, dtype=bool)

    selected_idx = np.where(select_mask)[0]
    if selected_idx.size == 0:
        raise ValueError("No variants selected from PLINK files")

    genotype_matrix = _read_plink_bed_subset(bed_path, n_samples, selected_idx)

    col_means = np.nanmean(genotype_matrix, axis=0)
    nan_r, nan_c = np.where(np.isnan(genotype_matrix))
    if nan_r.size:
        genotype_matrix[nan_r, nan_c] = col_means[nan_c]
    genotype_matrix = np.nan_to_num(genotype_matrix, nan=0.0)

    selected_bim = bim_df.iloc[selected_idx].reset_index(drop=True)
    rsid_out = selected_bim['RSID'].fillna('').tolist()
    chrpos_out = selected_bim['CHRPOS'].tolist()

    colnames = []
    for rsid_value, chrpos_value in zip(rsid_out, chrpos_out):
        if rsid_value:
            colnames.append(rsid_value)
        else:
            colnames.append(chrpos_value)

    if messages:
        warnings.warn("; ".join(messages))

    return {
        "matrix": genotype_matrix,
        "map": selected_bim[['CHR', 'RSID', 'POS', 'A1', 'A2']],
        "colnames": colnames,
        "chrpos": chrpos_out,
        "rsid": rsid_out
    }


# ----------------------------
# Variant / GWAS preparation
# ----------------------------

def load_variant_info(path: str):
    _require_pandas()
    df = pd.read_csv(path)
    if df.empty:
        return df
    lower = {c.lower(): c for c in df.columns}
    chr_col = next((lower[key] for key in ['chr', 'chrom', 'chromosome'] if key in lower), None)
    pos_col = next((lower[key] for key in ['pos_bp', 'pos', 'position', 'bp'] if key in lower), None)
    if chr_col is None or pos_col is None:
        raise ValueError('variant info must contain chr and position columns')
    df['chr'] = _canon_chr(df[chr_col])
    df['pos_bp'] = _canon_pos(df[pos_col])
    if 'rsid' not in df.columns:
        df['rsid'] = 'chr' + df['chr'].astype('Int64').astype(str) + ':' + df['pos_bp'].astype('Int64').astype(str)
    out = df[['rsid', 'chr', 'pos_bp']].dropna().drop_duplicates().sort_values(['chr', 'pos_bp'])
    return out.reset_index(drop=True)


def _detect_traits(df: "pd.DataFrame", trait_cols) -> List[str]:
    exclude = {'CHR', 'POS', 'CM', 'A1', 'A2', 'BP', 'N', 'SE', 'P', 'BETA'}
    numeric_cols = []
    for col in df.columns:
        if col.upper() in exclude:
            continue
        numeric = pd.to_numeric(df[col], errors='coerce')
        if numeric.notna().any():
            df[col] = numeric
            numeric_cols.append(col)
    if trait_cols:
        selected = [c for c in trait_cols if c in numeric_cols]
    else:
        selected = numeric_cols
    return selected[:2]


def auto_prepare_inputs(zscore_path: str, panel_prefix: str, trait_cols=None, verbose: bool = False) -> Dict[str, object]:
    _require_pandas()
    z_df = pd.read_csv(zscore_path, sep=None, engine='python')
    if 'CHR' not in z_df.columns or 'POS' not in z_df.columns:
        if 'RSID' in z_df.columns:
            split = z_df['RSID'].astype(str).str.split(':', expand=True)
            if split.shape[1] >= 2:
                z_df['CHR'] = split[0]
                z_df['POS'] = split[1]
    if 'CHR' not in z_df.columns or 'POS' not in z_df.columns:
        raise ValueError('zscore file must contain CHR and POS columns or rsid formatted as chr:pos')
    z_df['CHR'] = _canon_chr(z_df['CHR'])
    z_df['POS'] = _canon_pos(z_df['POS'])
    z_df = z_df.dropna(subset=['CHR', 'POS'])

    bim_path = f"{panel_prefix}.bim"
    bim_cols = ['CHR', 'RSID', 'CM', 'POS', 'A1', 'A2']
    bim_df = pd.read_csv(bim_path, sep=r'\s+', names=bim_cols)
    bim_df['CHR'] = _canon_chr(bim_df['CHR'])
    bim_df['POS'] = _canon_pos(bim_df['POS'])

    merged = pd.merge(bim_df, z_df, on=['CHR', 'POS'])
    if merged.empty:
        raise ValueError('No overlap between zscore and reference BIM on CHR+POS')

    traits = _detect_traits(merged, trait_cols)
    if not traits:
        raise ValueError('No numeric trait columns detected in zscore file')

    info_df = merged[['RSID', 'CHR', 'POS']].drop_duplicates().rename(columns={'POS': 'pos_bp'})
    info_df = info_df.sort_values(['CHR', 'pos_bp']).reset_index(drop=True)

    gwas_tables = []
    for trait in traits:
        trait_df = merged[['CHR', 'POS', trait]].rename(columns={trait: 'Z'})
        trait_df = trait_df.dropna(subset=['Z']).reset_index(drop=True)
        gwas_tables.append(trait_df)

    if verbose:
        print(f"auto_prepare_inputs: matched {info_df.shape[0]} variants; traits = {traits}")

    return {
        'variant_info': info_df,
        'traits': traits,
        'gwas_tables': gwas_tables
    }


def load_gwas_table(gwas_file: str) -> "pd.DataFrame":
    _require_pandas()
    df = pd.read_csv(gwas_file, sep=None, engine='python')
    if 'CHR' not in df.columns or 'POS' not in df.columns:
        raise ValueError('GWAS file must contain CHR and POS columns')
    df['CHR'] = _canon_chr(df['CHR'])
    df['POS'] = _canon_pos(df['POS'])
    z_col = next((c for c in df.columns if c.upper() == 'Z'), None)
    if z_col is None:
        raise ValueError('GWAS file must contain Z column')
    df['Z'] = pd.to_numeric(df[z_col], errors='coerce')
    df = df.dropna(subset=['CHR', 'POS', 'Z']).drop_duplicates(subset=['CHR', 'POS'])
    return df[['CHR', 'POS', 'Z']].reset_index(drop=True)


def load_multi_gwas_table(multi_file: str, zcols: str) -> Dict[str, object]:
    _require_pandas()
    df = pd.read_csv(multi_file, sep=None, engine='python')
    if 'CHR' not in df.columns or 'POS' not in df.columns:
        raise ValueError('multi_gwas file must contain CHR and POS columns')
    df['CHR'] = _canon_chr(df['CHR'])
    df['POS'] = _canon_pos(df['POS'])
    df = df.dropna(subset=['CHR', 'POS']).drop_duplicates(subset=['CHR', 'POS'])
    user_cols = [c.strip() for c in zcols.split(',') if c.strip()]
    for col in user_cols:
        if col not in df.columns:
            raise ValueError(f"z column '{col}' not found in multi_gwas file")
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return {
        'data': df,
        'zcols': user_cols
    }


def load_gene_catalog(coord_version: str, base_dir: str = "."):
    _require_pandas()
    candidate = coord_version.strip()
    if os.path.isfile(candidate):
        path = os.path.abspath(candidate)
    else:
        key = candidate.upper()
        if key not in GENE_FILES:
            raise FileNotFoundError(f"Unsupported gene catalog key or path: {coord_version}")
        filename = GENE_FILES[key]
        search_paths = [filename, os.path.join(base_dir, filename)]
        path = next((os.path.abspath(p) for p in search_paths if os.path.exists(p)), None)
        if path is None:
            raise FileNotFoundError(f"Cannot locate gene catalog file: {filename}")
    df = pd.read_csv(path, sep=None, engine='python')
    lower = {c.lower(): c for c in df.columns}
    rename_map = {
        'id': 'id',
        'gene_id': 'id',
        'symbol': 'id',
        'gene': 'id',
        'chr': 'chr',
        'chrom': 'chr',
        'chromosome': 'chr',
        'start': 'start',
        'start_bp': 'start',
        'position': 'start',
        'end': 'end',
        'stop': 'end',
        'stop_bp': 'end'
    }
    renamed = {}
    for key, target in rename_map.items():
        if key in lower and target not in renamed:
            renamed[target] = lower[key]
    df = df.rename(columns={orig: target for target, orig in renamed.items()})
    if 'id' not in df.columns:
        df['id'] = [f"gene_{i}" for i in range(1, len(df) + 1)]
    for col in ['chr', 'start', 'end']:
        if col not in df.columns:
            raise ValueError(f"gene catalog missing required column '{col}'")
    df['chr'] = df['chr'].astype(str).str.replace('^chr', '', regex=True, case=False)
    df['chr'] = pd.to_numeric(df['chr'], errors='coerce')
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['end'] = pd.to_numeric(df['end'], errors='coerce')
    df = df.dropna(subset=['id', 'chr', 'start', 'end'])
    df = df.sort_values(['chr', 'start', 'end']).reset_index(drop=True)
    return df[['id', 'chr', 'start', 'end']]


def annotate_nearest_gene(chromosomes, positions, gene_table):
    _require_pandas()
    import pandas as pd
    chrom_series = pd.to_numeric(pd.Series(chromosomes), errors='coerce')
    pos_series = pd.to_numeric(pd.Series(positions), errors='coerce')
    genes = gene_table.copy()
    if not isinstance(genes, pd.DataFrame):
        genes = pd.DataFrame(genes)
    required = {'id', 'chr', 'start', 'end'}
    if not required.issubset(set(genes.columns)):
        raise ValueError("gene_table must contain columns id, chr, start, end")
    genes = genes.copy()
    genes['chr'] = pd.to_numeric(genes['chr'], errors='coerce')
    genes['start'] = pd.to_numeric(genes['start'], errors='coerce')
    genes['end'] = pd.to_numeric(genes['end'], errors='coerce')
    genes = genes.dropna(subset=['chr', 'start', 'end'])
    genes = genes.sort_values(['chr', 'start', 'end'])

    results: List[Optional[str]] = [None] * len(pos_series)
    valid_mask = chrom_series.notna() & pos_series.notna()
    if not valid_mask.any():
        return results

    for chr_val in sorted(chrom_series[valid_mask].unique()):
        chr_mask = (chrom_series == chr_val) & valid_mask
        idx = np.where(chr_mask.to_numpy())[0]
        if idx.size == 0:
            continue
        genes_chr = genes[genes['chr'] == chr_val]
        if genes_chr.empty:
            continue
        starts = genes_chr['start'].to_numpy()
        ends = genes_chr['end'].to_numpy()
        ids = genes_chr['id'].astype(str).to_numpy()
        pos_vals = pos_series.iloc[idx].to_numpy()
        insert_idx = np.searchsorted(starts, pos_vals)
        prev_idx = np.clip(insert_idx - 1, 0, len(starts) - 1)
        next_idx = np.clip(insert_idx, 0, len(starts) - 1)
        dist_prev = np.where((pos_vals >= starts[prev_idx]) & (pos_vals <= ends[prev_idx]),
                             0,
                             np.minimum(np.abs(pos_vals - starts[prev_idx]), np.abs(pos_vals - ends[prev_idx])))
        dist_next = np.where((pos_vals >= starts[next_idx]) & (pos_vals <= ends[next_idx]),
                             0,
                             np.minimum(np.abs(pos_vals - starts[next_idx]), np.abs(pos_vals - ends[next_idx])))
        use_prev = dist_prev <= dist_next
        nearest = np.where(use_prev, ids[prev_idx], ids[next_idx])
        for arr_idx, gene_id in zip(idx, nearest):
            results[arr_idx] = gene_id
    return results


def align_genotype_to_gwas(geno_matrix,
                           geno_colnames,
                           geno_chrpos,
                           geno_rsid,
                           gwas_df,
                           prefer_rsid: bool = True):
    _require_pandas()
    import pandas as pd

    G = np.asarray(geno_matrix, dtype=np.float64)
    if G.ndim != 2:
        raise ValueError("geno_matrix must be 2D")

    colnames = list(geno_colnames) if geno_colnames is not None else None
    chrpos_attr = list(geno_chrpos) if geno_chrpos is not None else None
    rsid_attr = list(geno_rsid) if geno_rsid is not None else None

    gwas = gwas_df.copy()
    if 'chr' not in gwas.columns or 'pos_bp' not in gwas.columns:
        raise ValueError("gwas dataframe must contain chr and pos_bp columns")
    if 'rsid' not in gwas.columns:
        gwas['rsid'] = pd.NA
    gwas['chr'] = pd.to_numeric(gwas['chr'], errors='coerce')
    gwas['pos_bp'] = pd.to_numeric(gwas['pos_bp'], errors='coerce')
    gwas = gwas.dropna(subset=['chr', 'pos_bp'])
    gwas['chr'] = gwas['chr'].astype(int)
    gwas['pos_bp'] = gwas['pos_bp'].astype(int)
    gwas['chrpos'] = gwas['chr'].astype(str) + ':' + gwas['pos_bp'].astype(str)

    rsid_lookup = {}
    if rsid_attr is not None:
        for idx, val in enumerate(rsid_attr):
            if isinstance(val, str) and val not in rsid_lookup:
                rsid_lookup[val] = idx

    chrpos_lookup = {}
    if chrpos_attr is not None:
        for idx, val in enumerate(chrpos_attr):
            if isinstance(val, str) and val not in chrpos_lookup:
                chrpos_lookup[val] = idx

    colname_lookup = {}
    if colnames is not None:
        for idx, val in enumerate(colnames):
            if isinstance(val, str) and val not in colname_lookup:
                colname_lookup[val] = idx

    selected_cols: List[int] = []
    gwas_indices: List[int] = []
    variant_ids: List[str] = []
    positions: List[int] = []
    rsid_out: List[Optional[str]] = []

    for row_idx, row in enumerate(gwas.itertuples(index=False)):
        rsid_value = row.rsid if isinstance(row.rsid, str) else None
        chrpos_value = row.chrpos
        col_idx = None

        if prefer_rsid and rsid_value and rsid_value in rsid_lookup:
            col_idx = rsid_lookup[rsid_value]

        if col_idx is None and chrpos_value in chrpos_lookup:
            col_idx = chrpos_lookup[chrpos_value]

        if col_idx is None and colnames is not None and chrpos_value in colname_lookup:
            col_idx = colname_lookup[chrpos_value]

        if col_idx is None and not prefer_rsid and rsid_value and rsid_value in rsid_lookup:
            col_idx = rsid_lookup[rsid_value]

        if col_idx is None or col_idx in selected_cols:
            continue

        selected_cols.append(col_idx)
        gwas_indices.append(row_idx)
        variant_ids.append(chrpos_value)
        positions.append(int(row.pos_bp))
        rsid_out.append(rsid_value)

    if len(selected_cols) < 2:
        raise ValueError("Matched variants fewer than 2 between genotype matrix and GWAS data")

    G_selected = G[:, selected_cols]

    return {
        "matrix": G_selected,
        "variant_ids": variant_ids,
        "positions": positions,
        "gwas_indices": gwas_indices,
        "rsid": rsid_out,
        "matched": len(selected_cols)
    }
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


def qc_genotype_with_mask(G: "np.ndarray",
                          maf_threshold: float = 0.0,
                          mac_threshold: float = 25.0) -> Dict[str, object]:
    """Return filtered genotype matrix along with kept SNP indices."""
    if np is None:
        raise RuntimeError("numpy is required for accelerator")
    X = np.asarray(G, dtype=np.float64)
    X = np.where((X < 0) | (X > 2) | ~np.isfinite(X), 0.0, X)
    maf = X.mean(axis=0) / 2.0
    mac = X.sum(axis=0)
    var = X.var(axis=0)
    mask = (maf > maf_threshold) & (mac >= mac_threshold) & (var > 0.0) & np.isfinite(maf)
    indices = np.nonzero(mask)[0]
    filtered = X[:, indices]
    return {
        "matrix": filtered,
        "indices": indices
    }


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
