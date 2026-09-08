"""
Versión NO interactiva del control de calidad (QC) por muestra.

Reutiliza las mismas funciones y la misma lógica de máscara de
filtrado que el notebook original, pero no requiere intervención manual. 
Lee los umbrales ya decididos desde config/qc_thresholds.yaml y aplica la misma 
máscara de filtrado, con los mismos nombres de fichero de salida 

Este script es el que garantiza la reproducibilidad del TFM: cualquiera que
clone el repositorio y lo ejecute obtiene los mismos archivos filtrados que se
usaron en la memoria, sin depender de mover sliders correctamente.

"""

import gc
import glob
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


REPO_ROOT = Path(__file__).resolve().parents[1]
 
with open(REPO_ROOT / "config" / "params.yaml", "r", encoding="utf-8") as f:
    PARAMS = yaml.safe_load(f)
 
with open(REPO_ROOT / "config" / "qc_thresholds.yaml", "r", encoding="utf-8") as f:
    QC_CFG = yaml.safe_load(f)["qc"]
 
PROJECT_ROOT = os.environ.get(
    PARAMS["project_root_env_var"], PARAMS["project_root_default"]
)
 
RAW_H5AD_DIR = Path(PROJECT_ROOT) / PARAMS["paths"]["interim"]["00_raw_h5ad"]
FILT_DIR     = Path(PROJECT_ROOT) / PARAMS["paths"]["interim"]["01_qc_filtered"]
 
FIG_DIR = REPO_ROOT / PARAMS["paths"]["reports"]["figures"]["01_qc"]
REP_DIR = REPO_ROOT / PARAMS["paths"]["reports"]["tables"]["01_qc"]
 
for d in (FILT_DIR, FIG_DIR, REP_DIR):
    d.mkdir(parents=True, exist_ok=True)
 
MITO_ABSOLUTE_MAX = QC_CFG["mito_absolute_max"]
RUN_SCRUBLET = QC_CFG["run_scrublet"]
UMI_PREFILT_THRESHOLD_DEFAULT = QC_CFG["umi_prefilt_threshold_default"]
PER_SAMPLE = QC_CFG["per_sample"]
 
 
# =============================================================================
# Funciones auxiliares
# =============================================================================
 
def get_umi_counts(adata_raw):
    """Extrae el vector de UMIs totales por barcode (idéntico al original)."""
    if "counts" in adata_raw.layers:
        mat = adata_raw.layers["counts"]
        source = "layers['counts']"
    elif adata_raw.raw is not None:
        mat = adata_raw.raw.X
        source = "raw.X"
    else:
        mat = adata_raw.X
        source = "X"
 
    if sp.issparse(mat):
        umis = np.asarray(mat.sum(axis=1)).flatten()
    else:
        umis = np.asarray(mat.sum(axis=1)).flatten()
 
    return umis, source
 
 
def log10_ratio(adata):
    """Idéntico al original."""
    u = adata.obs["total_counts"].values
    g = adata.obs["n_genes_by_counts"].values
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.log10(u + 1) / np.log10(g + 1)
    return np.where(np.isfinite(r), r, 0.0)
 
 
def annotate_and_compute_qc(adata):
    """Idéntico al original (anotación MT/ribo + métricas + Scrublet)."""
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
 
    n_mt = adata.var["mt"].sum()
    n_ribo = adata.var["ribo"].sum()
 
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    )
    adata.obs["log10_ratio"] = log10_ratio(adata)
 
    if RUN_SCRUBLET:
        try:
            sc.pp.scrublet(adata, random_state=42)
        except Exception as e:
            print(f"  ⚠️  Scrublet falló: {e} — marcando todos como singletes.")
            adata.obs["predicted_doublet"] = False
            adata.obs["doublet_score"] = 0.0
    else:
        adata.obs["predicted_doublet"] = False
        adata.obs["doublet_score"] = 0.0
 
    return adata, n_mt, n_ribo
 
 
def plot_qc_panel(adata, sample_name, save_path=None):
    """Idéntico al panel original (violin + scatter + histogramas)."""
    obs = adata.obs
    C = dict(main="#2E86AB", mito="#E84855", ribo="#F4A261",
             ratio="#6A4C93", doub="#2A9D8F")
 
    fig = plt.figure(figsize=(20, 15))
    fig.suptitle(
        f"QC Metrics — {sample_name}\n"
        f"({adata.n_obs:,} células × {adata.n_vars:,} genes)",
        fontsize=13, fontweight="bold", y=0.99
    )
    gs = gridspec.GridSpec(3, 5, figure=fig, hspace=0.48, wspace=0.4)
 
    vmetrics = [
        ("n_genes_by_counts", "Genes / célula", C["main"]),
        ("total_counts", "UMIs / célula", C["main"]),
        ("pct_counts_mt", "% Mitocondrial", C["mito"]),
        ("pct_counts_ribo", "% Ribosómico", C["ribo"]),
        ("log10_ratio", "log₁₀(UMI)/log₁₀(g)", C["ratio"]),
    ]
    for i, (col, lbl, col_hex) in enumerate(vmetrics):
        ax = fig.add_subplot(gs[0, i])
        v = obs[col].dropna().values
        parts = ax.violinplot(v, showmedians=True, showextrema=True, widths=0.8)
        for pc in parts["bodies"]:
            pc.set_facecolor(col_hex); pc.set_alpha(0.7)
        parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(2)
        for k in ("cmaxes", "cmins", "cbars"):
            parts[k].set_color("#666")
        p5, p95 = np.percentile(v, [5, 95])
        ax.axhline(p5, color="orange", ls="--", lw=0.9, label=f"P5={p5:.0f}")
        ax.axhline(p95, color="red", ls="--", lw=0.9, label=f"P95={p95:.0f}")
        ax.set_title(lbl, fontsize=8, fontweight="bold")
        ax.set_xticks([]); ax.tick_params(axis="y", labelsize=7)
        ax.legend(fontsize=6, loc="upper right")
 
    ax2 = fig.add_subplot(gs[1, 0:2])
    sc_ = ax2.scatter(
        obs["total_counts"], obs["n_genes_by_counts"],
        c=obs["doublet_score"].values, cmap="RdYlGn_r",
        s=0.6, alpha=0.35, rasterized=True
    )
    plt.colorbar(sc_, ax=ax2, label="Doublet score", shrink=0.8)
    ax2.set_xlabel("UMIs totales", fontsize=8)
    ax2.set_ylabel("Genes detectados", fontsize=8)
    ax2.set_title("UMIs vs Genes\n(color = doublet score)", fontsize=9, fontweight="bold")
    ax2.tick_params(labelsize=7)
 
    ax3 = fig.add_subplot(gs[1, 2:4])
    ax3.scatter(obs["total_counts"], obs["pct_counts_mt"],
                c=C["mito"], s=0.6, alpha=0.35, rasterized=True)
    ax3.set_xlabel("UMIs totales", fontsize=8)
    ax3.set_ylabel("% Mitocondrial", fontsize=8)
    ax3.set_title("UMIs vs % Mitocondrial", fontsize=9, fontweight="bold")
    ax3.tick_params(labelsize=7)
 
    ax4 = fig.add_subplot(gs[1, 4])
    ax4.scatter(obs["pct_counts_mt"], obs["pct_counts_ribo"],
                c=C["doub"], s=0.6, alpha=0.35, rasterized=True)
    ax4.set_xlabel("% Mitocondrial", fontsize=8)
    ax4.set_ylabel("% Ribosómico", fontsize=8)
    ax4.set_title("% Mito vs % Ribo", fontsize=9, fontweight="bold")
    ax4.tick_params(labelsize=7)
 
    hcfg = [
        ("n_genes_by_counts", "Genes / célula", C["main"], None),
        ("total_counts", "UMIs / célula", C["main"], None),
        ("pct_counts_mt", "% Mitocondrial", C["mito"], (0, 60)),
        ("pct_counts_ribo", "% Ribosómico", C["ribo"], None),
        ("log10_ratio", "log₁₀ ratio", C["ratio"], (0.5, 1.5)),
    ]
    for i, (col, lbl, col_hex, xlim) in enumerate(hcfg):
        ax = fig.add_subplot(gs[2, i])
        v = obs[col].dropna().values
        ax.hist(v, bins=60, color=col_hex, alpha=0.8, edgecolor="none")
        for p, ls in zip([5, 25, 50, 75, 95], ["--", "--", "-", "--", "--"]):
            ax.axvline(np.percentile(v, p), color="black", ls=ls, lw=0.8, alpha=0.6)
        ax.set_xlabel(lbl, fontsize=8); ax.set_ylabel("Frecuencia", fontsize=8)
        ax.set_title(f"Distrib. {lbl}", fontsize=8, fontweight="bold")
        ax.tick_params(labelsize=7)
        if xlim:
            ax.set_xlim(xlim)
 
    stats = (
        f"mediana: genes={obs.n_genes_by_counts.median():.0f}  "
        f"UMIs={obs.total_counts.median():.0f}  "
        f"mito={obs.pct_counts_mt.median():.1f}%  "
        f"ribo={obs.pct_counts_ribo.median():.1f}%  "
        f"dobletes={obs.predicted_doublet.sum():,} "
        f"({100 * obs.predicted_doublet.mean():.1f}%)"
    )
    fig.text(0.5, 0.005, stats, ha="center", fontsize=8, style="italic", color="#444")
 
    if save_path:
        plt.savefig(save_path, dpi=110, bbox_inches="tight",
                    facecolor="white", edgecolor="none")
    plt.close(fig)
 
 
# =============================================================================
# Bucle principal aplica umbrales congelados
# =============================================================================
 
def process_sample(h5ad_path: Path, thresholds: dict) -> dict:
    sample_name = h5ad_path.stem
    print(f"\n{'=' * 60}\n  Muestra: {sample_name}\n{'=' * 60}")
 
    required = [
        "umi_prefilt_threshold", "min_genes", "max_genes", "min_umis",
        "max_umis", "max_mito_pct", "min_ribo_pct", "max_ribo_pct",
        "min_log10_ratio",
    ]
    missing = [k for k in required if thresholds.get(k) is None]
    if missing:
        raise ValueError(
            f"Faltan umbrales para '{sample_name}' en config/qc_thresholds.yaml: "
            f"{missing}. Rellénalos con los valores decididos en el notebook "
            f"interactivo antes de ejecutar este script."
        )
 
    adata_raw = sc.read_h5ad(h5ad_path)
    n_barcodes_total = adata_raw.n_obs
    umis_all, umi_source = get_umi_counts(adata_raw)
 
    thr_prefilt = thresholds["umi_prefilt_threshold"]
    mask_prefilt = umis_all >= thr_prefilt
    adata_qc = adata_raw[mask_prefilt].copy()
    n_cells_prefilt = adata_qc.n_obs
    del adata_raw
    gc.collect()
 
    adata_qc, n_mt, n_ribo = annotate_and_compute_qc(adata_qc)
    obs = adata_qc.obs
 
    fig_path_pre = FIG_DIR / f"{sample_name}_QC_prefiltro.png"
    plot_qc_panel(adata_qc, sample_name, save_path=fig_path_pre)
 
    mask = (
        (obs["n_genes_by_counts"] >= thresholds["min_genes"]) &
        (obs["n_genes_by_counts"] <= thresholds["max_genes"]) &
        (obs["total_counts"] >= thresholds["min_umis"]) &
        (obs["total_counts"] <= thresholds["max_umis"]) &
        (obs["pct_counts_mt"] <= thresholds["max_mito_pct"]) &
        (obs["pct_counts_ribo"] >= thresholds["min_ribo_pct"]) &
        (obs["pct_counts_ribo"] <= thresholds["max_ribo_pct"]) &
        (obs["log10_ratio"] >= thresholds["min_log10_ratio"])
    )
    if thresholds.get("remove_doublets", True):
        mask = mask & (~obs["predicted_doublet"].fillna(False).astype(bool))
 
    adata_filt = adata_qc[mask].copy()
    n_out = adata_filt.n_obs
    n_rem = adata_qc.n_obs - n_out
    pct_rem = 100 * n_rem / adata_qc.n_obs
 
    out_path = FILT_DIR / f"{sample_name}_filtered.h5ad"
    adata_filt.write_h5ad(out_path, compression="gzip")
 
    fig_path_post = FIG_DIR / f"{sample_name}_QC_postfiltro.png"
    plot_qc_panel(adata_filt, f"{sample_name} — POST filtro", save_path=fig_path_post)
 
    print(f"  Barcodes RAW        : {n_barcodes_total:,}")
    print(f"  Tras pre-filtro UMI  : {n_cells_prefilt:,} (≥ {thr_prefilt:,} UMIs)")
    print(f"  Tras QC definitivo   : {n_out:,} ({100 - pct_rem:.1f}% del pre-filtro)")
    print(f"  Guardado en          : {out_path}")
 
    stats_entry = {
        "sample": sample_name,
        "barcodes_raw": n_barcodes_total,
        "umi_source": umi_source,
        "umi_prefilt_thr": thr_prefilt,
        "cells_prefilt": n_cells_prefilt,
        "cells_filtered": n_out,
        "cells_removed_qc": n_rem,
        "pct_removed_qc": pct_rem,
        "min_genes": thresholds["min_genes"], "max_genes": thresholds["max_genes"],
        "min_umis": thresholds["min_umis"], "max_umis": thresholds["max_umis"],
        "max_mito_pct": thresholds["max_mito_pct"],
        "min_ribo_pct": thresholds["min_ribo_pct"], "max_ribo_pct": thresholds["max_ribo_pct"],
        "min_log10_ratio": thresholds["min_log10_ratio"],
        "remove_doublets": thresholds.get("remove_doublets", True),
        "genes_final": adata_filt.n_vars,
    }
 
    del adata_qc, adata_filt
    gc.collect()
    return stats_entry
 
 
def main():
    h5ad_files = sorted(glob.glob(str(RAW_H5AD_DIR / "*.h5ad")))
    if not h5ad_files:
        raise FileNotFoundError(f"No se encontraron archivos .h5ad en: {RAW_H5AD_DIR}")
 
    print(f"✅ Encontrados {len(h5ad_files)} archivos .h5ad en {RAW_H5AD_DIR}")
 
    all_stats = []
    for h5ad_path in h5ad_files:
        sample_name = Path(h5ad_path).stem
        thresholds = PER_SAMPLE.get(sample_name)
        if thresholds is None:
            raise KeyError(
                f"La muestra '{sample_name}' no existe en config/qc_thresholds.yaml "
                f"(sección qc.per_sample). Añádela antes de ejecutar el pipeline completo."
            )
        stats_entry = process_sample(Path(h5ad_path), thresholds)
        all_stats.append(stats_entry)
 
    df = pd.DataFrame(all_stats)
    csv_path = REP_DIR / "QC_summary_all_samples.csv"
    df.to_csv(csv_path, index=False, float_format="%.2f")
    print(f"\n✅ Resumen CSV guardado: {csv_path}")
 
    total_r = df["cells_prefilt"].sum()
    total_f = df["cells_filtered"].sum()
 
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("Resumen QC — Todas las muestras procesadas",
                 fontsize=13, fontweight="bold")
 
    x = range(len(df))
    xlbl = [s[:14] for s in df["sample"]]
 
    axes[0].bar(x, df["cells_prefilt"], label="Pre-filtro", color="#ADB5BD", alpha=0.9)
    axes[0].bar(x, df["cells_filtered"], label="Filtradas", color="#2E86AB", alpha=0.9)
    axes[0].set_xticks(x); axes[0].set_xticklabels(xlbl, rotation=45, ha="right", fontsize=7)
    axes[0].set_ylabel("Nº células"); axes[0].set_title("Células pre/post QC")
    axes[0].legend(fontsize=8)
 
    pct = df["pct_removed_qc"].astype(float)
    colors = ["#E84855" if p > 50 else "#F4A261" if p > 30 else "#2A9D8F" for p in pct]
    axes[1].bar(x, pct, color=colors, alpha=0.9)
    axes[1].axhline(30, color="orange", ls="--", lw=1, label="30%")
    axes[1].axhline(50, color="red", ls="--", lw=1, label="50% alerta")
    axes[1].set_xticks(x); axes[1].set_xticklabels(xlbl, rotation=45, ha="right", fontsize=7)
    axes[1].set_ylabel("% eliminado del pre-filtro")
    axes[1].set_title("% Células eliminadas en QC")
    axes[1].legend(fontsize=8)
 
    axes[2].pie(
        [total_f, total_r - total_f],
        labels=[f"Retenidas\n{total_f:,}", f"Eliminadas\n{total_r - total_f:,}"],
        colors=["#2E86AB", "#E84855"], autopct="%1.1f%%", startangle=90
    )
    axes[2].set_title(f"Total ({total_r:,} células candidatas)")
 
    plt.tight_layout()
    fig_path = REP_DIR / "QC_summary_global.png"
    plt.savefig(fig_path, dpi=120, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✅ Figura resumen guardada: {fig_path}")
 
    print(f"  Total pre-filtro : {total_r:,}")
    print(f"  Total filtradas  : {total_f:,} ({100 * total_f / total_r:.1f}% retenido)")
    print(f"  Muestras completadas: {len(df)}/{len(h5ad_files)}")
 
 
if __name__ == "__main__":
    main()
