# =========================================================================== #
# Figures C–E: CNA association analysis of treatment condition (Pre vs Post)
# within each NMF subtype. UMAP plots of neighborhood coefficients and
# lollipop charts of per-cluster average coefficients.
# =========================================================================== #

# ----------------------------- Import packages ------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
import cna
from scipy.stats import pearsonr
from matplotlib.patches import Patch

# ------------------------------ File paths ----------------------------------
input_h5ad        = "path/to/01merge.h5ad"                 # raw data
output_neighbors  = "path/to/01merge_neighbors.h5ad"       # processed with neighbors
output_umap_base  = "path/to/D_umap_case_coef"             # will append _{label}.pdf
output_lolli_base = "path/to/D_cluster_avg_coef_lollipop"  # will append _{label}.pdf

# ================= Part 1: Preprocess and save neighbors ====================
np.random.seed(0)

# Load data
d = ad.read_h5ad(input_h5ad)

# Extract sampleid and batch from index
d.obs['sampleid'] = d.obs.index.str.split('_').str[-1]
d.obs['batch']    = d.obs.index.str.split('_').str[-1]

# Keep rows with valid numeric identifiers
mask = (pd.to_numeric(d.obs['sampleid'], errors='coerce').notna()) & \
       (pd.to_numeric(d.obs['subtype'], errors='coerce').notna()) & \
       (pd.to_numeric(d.obs['seurat_clusters'], errors='coerce').notna())
d = d[mask].copy()

for col in ['sampleid', 'subtype', 'seurat_clusters', 'batch']:
    d.obs[col] = d.obs[col].astype('int64')

# Map treatment condition to integer
treatment_map = {'Unknow': 0, 'Pre': 1, 'Post': 2, 'Post2': 3}
d.obs['treatment_condition'] = d.obs['treatment_condition'].map(treatment_map).astype('int64')

# Compute PCA and neighbors
if 'X_pca' not in d.obsm:
    sc.pp.pca(d)
sc.pp.neighbors(d, n_neighbors=40)

# Remove .raw to avoid write issues
if d.raw is not None:
    d.raw = None
if '_index' in d.var.columns:
    d.var = d.var.drop(columns='_index')

d.write_h5ad(output_neighbors)

# ================= Part 2: Subtype-specific CNA analysis ====================
d_all = ad.read_h5ad(output_neighbors)
# Keep only Pre (1) and Post (2)
d_all = d_all[~d_all.obs['treatment_condition'].isin([0, 3]), :].copy()

subtype_map = {1: 'EC1', 2: 'EC2', 3: 'EC3'}

for sub, label in subtype_map.items():
    # Subset to current subtype
    exclude = [s for s in [1, 2, 3] if s != sub]
    d = d_all[~d_all.obs['subtype'].isin(exclude), :].copy()

    # Sample-level metadata
    samplem = cna.ut.obs_to_sample(
        d,
        ['subtype', 'seurat_clusters', 'treatment_condition', 'batch'],
        'sampleid'
    )
    samplem['case_status'] = (samplem['treatment_condition'] == 2).astype(int)

    # Run CNA association
    res = cna.tl.association(
        d,
        samplem['case_status'],
        'sampleid',
        batches=None,
        key_added='case_coef',
        return_full=True
    )

    # Compute raw p-values (fallback to manual Pearson if not provided)
    if hasattr(res, 'pvals'):
        raw_pvals = res.pvals
    elif hasattr(res, 'pvalues'):
        raw_pvals = res.pvalues
    else:
        raw_pvals = None
        if hasattr(res, 'namresid'):
            X = res.namresid.values
            sample_ids = res.namresid.index
            y_vals = samplem.loc[sample_ids, 'case_status'].values.astype(float)
            n_nei = X.shape[1]
            raw_pvals = np.empty(n_nei)
            for i in range(n_nei):
                _, p_val = pearsonr(X[:, i], y_vals)
                raw_pvals[i] = p_val

    if raw_pvals is not None:
        d.obs['case_coef_pval'] = np.nan
        d.obs.loc[res.kept, 'case_coef_pval'] = raw_pvals

    # ---------- UMAP plot (all neighborhoods) ----------
    umap_coords = d.obsm['X_umap']
    coef_vals = d.obs['case_coef'].values

    vmin, vmax = -0.2, 0.2
    cmap = plt.cm.RdBu_r
    norm = plt.matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap.set_over('red')
    cmap.set_under('blue')

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
               c=coef_vals, cmap=cmap, norm=norm, s=0.5,
               rasterized=True, zorder=0)
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, extend='both')
    cbar.set_label('Association coefficient (case_coef)')
    ax.set_title(f'UMAP: case_coef (all neighborhoods) - {label}')
    plt.savefig(f'{output_umap_base}_{label}.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    # ---------- Per-cluster average coefficient lollipop ----------
    df = d.obs[['seurat_clusters', 'case_coef']].copy()
    df['seurat_clusters'] = df['seurat_clusters'].astype(str)

    cluster_avg = df.groupby('seurat_clusters')['case_coef'].mean().reset_index()
    cluster_avg.columns = ['cluster', 'avg_coef']
    cluster_avg = cluster_avg[np.abs(cluster_avg['avg_coef']) > 0.1]
    cluster_avg = cluster_avg.sort_values('avg_coef', ascending=False)

    # Scale factor for EC3 (smaller points)
    scale = 0.3 if label == 'EC3' else 1.0

    x = np.arange(len(cluster_avg))
    y = cluster_avg['avg_coef'].values
    labels = cluster_avg['cluster'].values

    size = np.clip((np.abs(y) ** 2) * 3000 * scale, 10 * scale, 400 * scale)
    point_colors = ['#d62728' if val > 0 else '#1f77b4' for val in y]

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.scatter(x, y, s=size, c=point_colors, alpha=0.8, edgecolors='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, ha='right', fontsize=10)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.set_xlabel('Seurat Clusters')
    ax.set_ylabel('Average case_coef')
    ax.set_title(f'Average Association Coefficient per Cluster (|avg_coef| > 0.1) - {label}\n'
                 '(Point size represents |average coefficient|)')

    # Color legend
    legend_color = [Patch(facecolor='#d62728', edgecolor='black', label='Positive coef'),
                    Patch(facecolor='#1f77b4', edgecolor='black', label='Negative coef')]
    leg1 = ax.legend(handles=legend_color, loc='upper left', bbox_to_anchor=(1.02, 1.0),
                     title='Coefficient sign', frameon=True)
    ax.add_artist(leg1)

    # Size legend
    fixed_abs_vals = [0.1, 0.2, 0.3, 0.4]
    abs_y = np.abs(y)
    if len(abs_y) > 0:
        min_abs, max_abs = abs_y.min(), abs_y.max()
        size_vals_to_show = [v for v in fixed_abs_vals if v >= min_abs and v <= max_abs]
    else:
        size_vals_to_show = []

    size_handles = []
    for val in size_vals_to_show:
        s = np.clip(val * 500 * scale, 20 * scale, 200 * scale)
        handle = ax.scatter([], [], s=s, c='gray', alpha=0.8,
                            edgecolors='black', linewidth=0.5,
                            label=f'{val:.1f}')
        size_handles.append(handle)
    ax.legend(handles=size_handles, loc='lower left', bbox_to_anchor=(1.02, 0.0),
              title='|avg_coef|', frameon=True)

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(f'{output_lolli_base}_{label}.pdf', dpi=300, bbox_inches='tight')
    plt.close()