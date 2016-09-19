#!/usr/bin/env python


"""
Usage:

sbatch -p develop --mem 80000 -c 4 -J qc_plots \
-o /scratch/lab_bock/shared/projects/crop-seq/results_pipeline/qc_plots.log \
~/qc_plots.sh

"""

import matplotlib
matplotlib.use('pdf')
import argparse
import os
import pandas as pd
import pysam
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


root_dir = "/scratch/lab_bock/shared/projects/crop-seq"
results_dir = os.path.join(root_dir, "results")
sample_annotation = pd.read_csv(os.path.join(root_dir, "metadata/annotation.csv"))

# get guide annotation
guide_annotation = os.path.join(root_dir, "metadata/guide_annotation.csv")
guide_annotation = pd.read_csv(guide_annotation)

reads = pd.read_csv(os.path.join(results_dir, "guide_cell_quantification.all.csv"))
scores = pd.read_csv(os.path.join(results_dir, "guide_cell_scores.all.csv"))
assignment = pd.read_csv(os.path.join(results_dir, "guide_cell_assignment.all.csv"))
# group gRNAs per gene or control group
assignment = assignment.set_index("assignment")
assignment["group"] = pd.np.nan
assignment.loc[assignment.index.str.contains("CTRL"), "group"] = "CTRL"
assignment.loc[assignment.index.str.contains("Essential"), "group"] = "Essential"
assignment = assignment.reset_index()
tr = assignment[assignment["group"].isnull()]['assignment'].str.extract(".*_(.*)_.*")
assignment.loc[tr.index, "group"] = tr

# group by experiment (without stimulation condition)
assignment["experiment_group"] = assignment["experiment"].str.extract("^(.*)_.*$")


# Stats numbers
all_annot = pd.read_csv(os.path.join(root_dir, "metadata/annotation.all.csv"))
all_annot["merged_sample_name"] = "CROP-seq_" + all_annot["merged_sample_name"]

stats = pd.DataFrame()
for sample_name in sample_annotation[~sample_annotation["grna_library"].isnull()]["sample_name"].unique():
    # Processing stats
    if "WNT" in sample_name:
        s = pd.read_csv(os.path.join(root_dir, "results_pipeline", sample_name, "stats.tsv"), header=None, names=["metric", "value", "pipeline"], sep="\t")

        # Remove strings
        s = s[s["metric"] != "Success"]
        s = s[s["metric"] != "Time"]
        s = s[~s["metric"].str.contains(" on")]
        # Colapse some
        s['metric'] = s['metric'].str.replace("\(.*", "")
        # Values to numerical
        s['value'] = s['value'].str.strip().str.strip("%").str.strip("\|").astype(float)
    else:
        ss = pd.DataFrame()
        for name in all_annot[all_annot["merged_sample_name"] == sample_name]["sample_name"]:
            print(sample_name, name)
            ds = pd.read_csv(os.path.join(root_dir, "results_pipeline", name, "stats.tsv"), header=None, names=["metric", "value", "pipeline"], sep="\t")

            # Remove strings
            ds = ds[ds["metric"] != "Success"]
            ds = ds[ds["metric"] != "Time"]
            ds = ds[~ds["metric"].str.contains(" on")]
            # Colapse some
            ds['metric'] = ds['metric'].str.replace("\(.*", "")
            # Values to numerical
            ds['value'] = ds['value'].str.strip().str.strip("%").str.strip("\|").astype(float)
            ss = ss.append(ds)

        # Sum up metrics
        s = ss.groupby("metric")['value'].apply(sum)
        # or average if percentages
        s.loc[s.index[s.index.str.contains("%|per ")], ] = s[s.index[s.index.str.contains("%|per |verage")]] / 4.
        s = s.reset_index()
        s.columns = ["metric", "value"]
        s["pipeline"] = "dropseq"

    # Number of assigned cells
    c = pd.read_csv(os.path.join(root_dir, "results_pipeline", sample_name, "quantification", "guide_cell_assignment.csv"))['cell'].drop_duplicates().shape[0]
    s = s.append({"metric": "number_assigned_cells", "value": c, "pipeline": "dropseq"}, ignore_index=True)

    s["sample_name"] = sample_name
    stats = stats.append(s)

stats.to_csv(os.path.join(results_dir, "pipeline_stats.csv"), index=False)
stats = pd.read_csv(os.path.join(results_dir, "pipeline_stats.csv"))

stats_pivot = pd.pivot_table(stats, index="metric", columns="sample_name", values="value")
stats_pivot = pd.read_csv(os.path.join(results_dir, "pipeline_stats.pivot.csv"), index_col=0)

# Some cell-centric stats now:
for n_genes in [500]:
    for sample_name in sample_annotation[~sample_annotation["grna_library"].isnull()]["sample_name"].unique()[:3]:
        exp = pd.read_csv(os.path.join(root_dir, "results_pipeline", sample_name, "digital_expression.{}genes.tsv".format(n_genes)), index_col=0, sep="\t")
        print(sample_name, exp.shape)

        # metrics:
        # cells_per_gene = exp.apply(lambda x: (x > 0).sum(), axis=1)
        # reads per cell
        reads_per_cell = exp.sum(axis=0)
        stats_pivot.loc[" ".join([str(n_genes), "unique used reads"]), sample_name] = reads_per_cell.sum()

        # genes per cell
        genes_per_cell = exp.apply(lambda x: (x > 0).sum(), axis=0)
        # % mitochondrial
        mito_per_cell = (exp.ix[exp.index[exp.index.str.contains("MT-")]].sum(axis=0) / reads_per_cell) * 100
        # % ribosomal proteins
        ribo_per_cell = (exp.ix[exp.index[exp.index.str.contains("RP")]].sum(axis=0) / reads_per_cell) * 100

        # for each, add to stats table distribution values (mean, 5th, 25th, 50th, 75th, 95th percentiles)
        for d in ["reads_per_cell", "genes_per_cell", "mito_per_cell", "ribo_per_cell"]:
            for metric in reads_per_cell.describe().index:
                stats_pivot.loc[" ".join([str(n_genes), d, metric]), sample_name] = eval(d).describe().ix[metric]

        # plot distributions of each
        fig, axis = plt.subplots(2, 2, figsize=(12, 8), sharex=False, sharey=False)
        axis = axis.flatten()
        for i, d in enumerate(["reads_per_cell", "genes_per_cell", "mito_per_cell", "ribo_per_cell"]):
            sns.distplot(eval(d), ax=axis[i])
            axis[i].set_title(d)
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "cell_stats.{}.{}.distplots.svg".format(sample_name, n_genes)), bbox_inches="tight")

        # plot pairwise distributions of each
        combined_stats = pd.concat([reads_per_cell, genes_per_cell, mito_per_cell, ribo_per_cell], 1)
        combined_stats.columns = ["reads_per_cell", "genes_per_cell", "mito_per_cell", "ribo_per_cell"]
        combined_stats.to_csv(os.path.join(results_dir, "cell_stats.{}.{}.csv".format(sample_name, n_genes)), index=True)
        g = sns.pairplot(combined_stats, diag_kws={"bins": 100}, plot_kws={"s": 10, "alpha": 0.5})
        g.fig.savefig(os.path.join(results_dir, "cell_stats.{}.{}.pairwise_dists.svg".format(sample_name, n_genes)), bbox_inches="tight")

        #

        # UMI duplication stats
        umi = pd.read_csv(os.path.join(root_dir, "results_pipeline", sample_name, "cell_umi_barcodes.{}genes.tsv".format(n_genes)), sep="\t")

        # percent unique UMIs and UMI duplication
        stats_pivot.loc[" ".join([str(n_genes), "percent unique UMIs"]), sample_name] = ((umi["Num_Obs"] == 1).sum() / umi.shape[0]) * 100
        for n in range(2, 9):
            stats_pivot.loc[" ".join([str(n_genes), "percent UMIs in {} copies".format(n)]), sample_name] = ((umi["Num_Obs"] == n).sum() / umi.shape[0]) * 100
        stats_pivot.loc[" ".join([str(n_genes), "percent UMIs in >8 copies"]), sample_name] = ((umi["Num_Obs"] > 8).sum() / umi.shape[0]) * 100
        # total number of used reads (unique reads across all exons, all genes, all cells)
        stats_pivot.loc[" ".join([str(n_genes), "total unique used reads"]), sample_name] = reads_per_cell.sum()
        stats_pivot.to_csv(os.path.join(results_dir, "pipeline_stats.pivot.csv"), index=True)

        # plot distibution
        fig, axis = plt.subplots(1, figsize=(12, 8))
        sns.distplot(np.log10(umi["Num_Obs"]), kde=False, ax=axis)
        axis.set_yscale('log')
        axis.set_ylabel("Frequency")
        axis.set_xlabel("UMI observations (log10)")
        fig.savefig(os.path.join(results_dir, "umi_stats.{}.{}.distplot.svg".format(sample_name, n_genes)), bbox_inches="tight")

        #

        # gRNA assignemnt
        # number of cells in transcriptome matrix that are assigned

        # total
        exp_assigments = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in exp.columns])
        stats_pivot.loc[" ".join([str(n_genes), "total assigned cells"]), sample_name] = sum([type(x) != str for x in exp_assigments])
        stats_pivot.loc[" ".join([str(n_genes), "total unassigned cells"]), sample_name] = sum([type(x) == str for x in exp_assigments])

        # correctly assigned (assigned to gRNA of proper library)
        if "WNT" in sample_name:
            cor_assignment = assignment[assignment["assignment"].str.contains("Wnt")]
        elif "TCR" in sample_name:
            cor_assignment = assignment[assignment["assignment"].str.contains("Tcr")]

        exp_assigments = pd.Series([cor_assignment[cor_assignment["cell"] == y]["group"].squeeze() for y in exp.columns])
        stats_pivot.loc[" ".join([str(n_genes), "correctly assigned cells"]), sample_name] = sum([type(x) != str for x in exp_assigments])
        stats_pivot.loc[" ".join([str(n_genes), "correctly unassigned cells"]), sample_name] = sum([type(x) == str for x in exp_assigments])

        # plot number of genes assigned depending on the minimum number of genes required
        total = list()
        assigned = list()
        for j in range(1, 10000, 10):
            print(j)
            pos = genes_per_cell > j
            d = genes_per_cell[pos]
            try:
                original_names = pd.Series(d.index.str.split("_")).apply(lambda x: x[0])
                assigned.append(original_names.isin(scores["cell"].dropna()).sum())
                total.append(pos.sum())
            except AttributeError:
                break

        fig, axis = plt.subplots(2, sharex=False, sharey=False, figsize=(8, 8))
        axis[0].plot((range(1, len(total) * 10, 10)), (total), lw=4, color="b", label="All cells")
        axis[0].plot((range(1, len(total) * 10, 10)), (assigned), lw=4, color="r", label="Cells with gRNA assigned")
        axis[1].plot(np.log2(range(1, len(total) * 10, 10)), np.log2(total), lw=4, color="b", label="All cells")
        axis[1].plot(np.log2(range(1, len(total) * 10, 10)), np.log2(assigned), lw=4, color="r", label="Cells with gRNA assigned")
        axis[0].set_xlabel("Genes detected")
        axis[1].set_xlabel("Genes detected (log2)")
        axis[0].set_ylabel("Cells")
        axis[1].set_ylabel("Cells (log2)")
        plt.legend()
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "barcodes_per_cell.{}.{}.identified.varying_genes_detected.svg".format(sample_name, n_genes)), bbox_inches="tight")

        # cumulative sum
        d = pd.DataFrame([total, assigned], index=['total', 'assigned'], columns=range(1, len(total) * 10, 10)).T
        cumsum = d.cumsum() / d.sum()

        fig, axis = plt.subplots(2, sharex=False, sharey=False, figsize=(8, 8))
        axis[0].plot(cumsum.index, cumsum['total'], lw=4, color="b", label="All cells")
        axis[0].plot(cumsum.index, cumsum['assigned'], lw=4, color="r", label="Cells with gRNA assigned")
        axis[1].plot(np.log2(cumsum.index), cumsum['total'], lw=4, color="b", label="All cells")
        axis[1].plot(np.log2(cumsum.index), cumsum['assigned'], lw=4, color="r", label="Cells with gRNA assigned")
        axis[0].set_xlabel("Genes detected")
        axis[1].set_xlabel("Genes detected (log2)")
        axis[0].set_ylabel("Cells (fraction of total)")
        axis[1].set_ylabel("Cells (fraction of total)")
        plt.legend()
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "barcodes_per_cell.{}.{}.identified.varying_genes_detected.cumulative.svg".format(sample_name, n_genes)), bbox_inches="tight")

        # Plot reads vs genes covered for each cell
        reads_per_cell = exp.sum(0)

        fig, axis = plt.subplots(2, sharex=False, sharey=False, figsize=(8, 8))
        axis[0].scatter(reads_per_cell, genes_per_cell, s=4)
        axis[1].scatter(np.log2(reads_per_cell), np.log2(genes_per_cell), s=4)
        axis[0].set_xlabel("Reads")
        axis[1].set_xlabel("Reads (log2)")
        axis[0].set_ylabel("Genes")
        axis[1].set_ylabel("Genes (log2)")
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "reads_vs_genes.{}.{}.svg".format(sample_name, n_genes)), bbox_inches="tight")


metrics = [
    "TagBamWithReadSequenceExtended_cell in total ",
    "TagBamWithReadSequenceExtended_cell paired in sequencing",
    "PolyATrimmer in total ",
    "TrimStartingSequence in total ",
    "FilterBAM in total ",
    "STAR Number of input reads |",
    "STAR Average input read length |",
    "STAR Average mapped length |",
    "STAR Mismatch rate per base, % |",
    "STAR Number of splices: Total |",
    "STAR % of reads unmapped: too many mismatches |",
    "STAR % of reads unmapped: too short |",
    "STAR % of reads unmapped: other |",
    "STAR Number of reads mapped to multiple loci |",
    "STAR % of reads mapped to multiple loci |",
    "STAR Number of reads mapped to too many loci |",
    "STAR % of reads mapped to too many loci |",
    "STAR Uniquely mapped reads number |",
    "STAR Uniquely mapped reads % |",
    "BAMTagHistogram in total ",
    "BAMTagHistogram mapped ",
    "SortSam mapped ",
    "MergeBamAlignment in total ",
    "TagReadWithGeneExon mapped ",
    "500 percent unique UMIs",
    "500 percent UMIs in 2 copies",
    "500 percent UMIs in 3 copies",
    "500 percent UMIs in 4 copies",
    "500 percent UMIs in 5 copies",
    "500 percent UMIs in 6 copies",
    "500 percent UMIs in 7 copies",
    "500 percent UMIs in 8 copies",
    "500 percent UMIs in >8 copies",
    "500 total unique used reads",
    "DigitalExpression_500genes number_cells",
    "DigitalExpression_500genes number_genes",
    "500 total assigned cells",
    "500 total unassigned cells",
    "500 correctly assigned cells",
    "500 correctly unassigned cells",
    "number_assigned_cells",
    "500 reads_per_cell 25%",
    "500 reads_per_cell 50%",
    "500 reads_per_cell mean",
    "500 reads_per_cell 75%",
    "500 genes_per_cell 25%",
    "500 genes_per_cell 50%",
    "500 genes_per_cell mean",
    "500 genes_per_cell 75%",
    "500 mito_per_cell 25%",
    "500 mito_per_cell 50%",
    "500 mito_per_cell mean",
    "500 mito_per_cell 75%",
    "500 ribo_per_cell 25%",
    "500 ribo_per_cell 50%",
    "500 ribo_per_cell mean",
    "500 ribo_per_cell 75%"]

stats_pivot = stats_pivot.ix[metrics]
# stats_pivot.index = stats_pivot.index.str.extract("(.*) |") #.str.strip()
stats_pivot.to_csv(os.path.join(results_dir, "pipeline_stats.pivot.selected.csv"), index=True)

# plot selected metrics
fig, axis = plt.subplots(2, figsize=(12, 16))
sns.heatmap(stats_pivot, ax=axis[0], annot=True, fmt='.1f')
sns.heatmap(stats_pivot.apply(lambda x: (x - x.mean()) / x.std(), axis=1), ax=axis[1], annot=False)
plt.setp(axis[0].xaxis.get_majorticklabels(), visible=False)
plt.setp(axis[1].xaxis.get_majorticklabels(), rotation=90)
plt.setp(axis[0].yaxis.get_majorticklabels(), rotation=0)
plt.setp(axis[1].yaxis.get_majorticklabels(), rotation=0)
fig.savefig(os.path.join(results_dir, "pipeline_stats.pivot.selected.svg"), bbox_inches="tight")
