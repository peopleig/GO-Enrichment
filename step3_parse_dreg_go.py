import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from gprofiler import GProfiler

DREG_FILE = "data/nrf1_dreg_hitss.dreg"
HITS_CSV = "results/nrf1_hit_genes.csv"
GO_CSV = "results/go_enrichment_results.csv"
PLOT_BP = "results/go_barplot_BP.png"
PLOT_MF = "results/go_barplot_MF.png"
PLOT_CC = "results/go_barplot_CC.png"
PLOT_DOTPLOT = "results/go_dotplot.png"

import os
os.makedirs("results", exist_ok=True)

# --- Parse dreg output to get gene names ---
genes = set()
with open(DREG_FILE) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # dreg excel-format: SeqName \t Start \t End \t Score \t Strand
        # SeqName looks like: chr1@11869-11870|DDX11L2::chr1:11369-11870(+)
        # We parse the gene name from the | delimiter
        parts = line.split("\t")
        if len(parts) < 1:
            continue
        seq_name = parts[0]
        match = re.search(r"\|([^|:]+?)(?:::|$)", seq_name)
        if match:
            genes.add(match.group(1).strip())

genes = sorted(genes)
print(f"Genes with NRF1 motif in promoter: {len(genes)}")

pd.DataFrame({"gene": genes}).to_csv(HITS_CSV, index=False)
print(f"Saved gene list to {HITS_CSV}")

# --- GO enrichment via g:Profiler ---
gp = GProfiler(return_dataframe=True)
go_results = gp.profile(
    organism="hsapiens",
    query=genes,
    sources=["GO:BP", "GO:MF", "GO:CC"],
    significance_threshold_method="fdr",
    user_threshold=0.05,
    no_evidences=False,
)

if go_results.empty:
    print("No significant GO terms found.")
    sys.exit(0)

go_results = go_results.sort_values("p_value")
go_results.to_csv(GO_CSV, index=False)
print(f"Saved GO results to {GO_CSV} ({len(go_results)} terms)")

# --- Plotting helper ---
def make_barplot(df, source_label, outfile, top_n=20, color_base="#4a7fb5"):
    subset = df[df["source"] == source_label].head(top_n).copy()
    if subset.empty:
        print(f"No results for {source_label}, skipping barplot.")
        return
    subset = subset.sort_values("p_value", ascending=False)
    subset["-log10(p)"] = -np.log10(subset["p_value"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(8, max(4, len(subset) * 0.35)))
    colors = cm.Blues(np.linspace(0.4, 0.85, len(subset)))
    bars = ax.barh(subset["name"], subset["-log10(p)"], color=colors)
    ax.set_xlabel("-log10(adjusted p-value)")
    ax.set_title(f"GO Enrichment — {source_label} (top {len(subset)})")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {outfile}")

def make_dotplot(df, outfile, top_n=30):
    subset = df.head(top_n).copy()
    if subset.empty:
        return
    subset = subset.sort_values("p_value", ascending=False)
    subset["-log10(p)"] = -np.log10(subset["p_value"].clip(lower=1e-300))
    subset["gene_ratio"] = subset["intersection_size"] / subset["query_size"]

    source_colors = {"GO:BP": "#e07b39", "GO:MF": "#4a7fb5", "GO:CC": "#5aa45a"}
    point_colors = subset["source"].map(source_colors).fillna("gray")

    fig, ax = plt.subplots(figsize=(8, max(5, len(subset) * 0.35)))
    sc = ax.scatter(
        subset["gene_ratio"],
        subset["name"],
        c=subset["-log10(p)"],
        s=subset["intersection_size"] * 3,
        cmap="YlOrRd",
        edgecolors="gray",
        linewidths=0.4,
        alpha=0.85,
    )
    cbar = plt.colorbar(sc, ax=ax, shrink=0.5)
    cbar.set_label("-log10(p-value)")
    ax.set_xlabel("Gene ratio")
    ax.set_title(f"GO Enrichment dotplot (top {len(subset)} terms)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {outfile}")

make_barplot(go_results, "GO:BP", PLOT_BP)
make_barplot(go_results, "GO:MF", PLOT_MF)
make_barplot(go_results, "GO:CC", PLOT_CC)
make_dotplot(go_results, PLOT_DOTPLOT)

print("\nAll done. Outputs in results/")