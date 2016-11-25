#!/usr/bin/env python

import os
import pandas as pd
import pysam
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from looper.models import Project


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def get_reads_in_construct(bam, guide_annotation):
    def overlap_1d(min1, max1, min2, max2):
        return max(0, min(max1, max2) - max(min1, min2))

    bam_handle = pysam.AlignmentFile(bam)

    reads = pd.DataFrame()
    # for each "chromosome" (each guideRNA)
    for chrom in guide_annotation["target"].unique():
        print(chrom)

        if chrom == "Cas9_blast":
            continue

        # get position of alignment
        guide_seq = guide_annotation[guide_annotation["target"] == chrom]['sequence'].squeeze()
        chrom_size = len(prj.config['crop-seq']['u6'] + guide_seq + prj.config['crop-seq']['rest'])
        guide_start_pos = len(prj.config['crop-seq']['u6']) + 1
        guide_end_pos = chrom_size - len(prj.config['crop-seq']['rest'])

        # for each read
        for aln in bam_handle.fetch(reference=chrom + "_chrom"):
            # skip reads
            if (
                aln.is_qcfail or  # failed quality (never happens, but for the future)
                aln.is_secondary or
                np.mean(aln.query_alignment_qualities) < 10 or  # low mapping Q (never happens, but for the future)
                "--" in aln.get_reference_sequence()  # reads with two+ gaps
            ):
                continue

            # get cell index
            cell = dict(aln.get_tags())['XC']

            # get molecule index
            molecule = dict(aln.get_tags())['XM']

            # determine distance to end of gRNA sequence
            distance = aln.reference_start - guide_end_pos

            # determine numbner of overlaping bases
            overlap = overlap_1d(aln.reference_start, aln.reference_end, guide_start_pos, guide_end_pos)

            # determine if inside gRNA
            inside = True if overlap > 0 else False
            # make sure strand is correct
            if aln.is_reverse:
                strand_agreeement = False
            else:
                strand_agreeement = True
            # get alignement quality
            # aln.mapapping_quality
            mapping_quality = np.mean(aln.query_alignment_qualities)

            reads = reads.append(pd.Series([chrom, cell, molecule, distance, overlap, inside, mapping_quality, strand_agreeement]), ignore_index=True)
    reads.columns = ["chrom", "cell", "molecule", "distance", "overlap", "inside", "mapping_quality", "strand_agreeement"]
    return reads


def get_reads_in_Cas9_construct(bam):
    def overlap_1d(min1, max1, min2, max2):
        return max(0, min(max1, max2) - max(min1, min2))

    chrom = "Cas9_blast"

    bam_handle = pysam.AlignmentFile(bam)

    reads = pd.DataFrame()

    sequence = "".join([
        prj.config['crop-seq']['cas9'],
        prj.config['crop-seq']['nls'],
        prj.config['crop-seq']['flag'],
        prj.config['crop-seq']['p2a'],
        prj.config['crop-seq']['blast'],
        prj.config['crop-seq']['space'],
        prj.config['crop-seq']['virus_ltr']])

    # get position of alignment
    chrom_size = len(sequence)
    guide_start_pos = 0
    guide_end_pos = chrom_size - len(prj.config['crop-seq']['cas9'])

    # for each read
    for aln in bam_handle.fetch(reference=chrom + "_chrom"):
        # skip reads
        if (
            aln.is_qcfail or  # failed quality (never happens, but for the future)
            aln.is_secondary or
            np.mean(aln.query_alignment_qualities) < 10 or  # low mapping Q (never happens, but for the future)
            "--" in aln.get_reference_sequence()  # reads with two+ gaps
        ):
            continue

        # determine distance to start of Cas9 construct
        distance = guide_start_pos - aln.reference_start

        if distance < 0:
            continue

        # determine numbner of overlaping bases
        overlap = overlap_1d(aln.reference_start, aln.reference_end, guide_start_pos, guide_end_pos)

        # determine if overlaps Cas9 construct
        inside = True if overlap > 0 else False
        # make sure strand is correct
        if aln.is_reverse:
            strand_agreeement = False
        else:
            strand_agreeement = True
        # get alignement quality
        # aln.mapapping_quality
        mapping_quality = np.mean(aln.query_alignment_qualities)

        # get cell index
        cell = dict(aln.get_tags())['XC']

        # get molecule index
        molecule = dict(aln.get_tags())['XM']

        reads = reads.append(pd.Series([chrom, cell, molecule, distance, overlap, inside, mapping_quality, strand_agreeement]), ignore_index=True)
    reads.columns = ["chrom", "cell", "molecule", "distance", "overlap", "inside", "mapping_quality", "strand_agreeement"]
    return reads


def plot_reads_in_constructs(reads):
    from collections import Counter
    # Inspect
    fig, axis = plt.subplots(2, sharex=True)
    # number of barcode reads per cell
    sns.distplot(np.log2(1 + reads.groupby(["cell"]).apply(len)), ax=axis[0], kde=False)
    axis[0].set_xlabel("Reads (log2(1 + x))")
    # number of unique barcode reads per cell
    sns.distplot(np.log2(1 + reads.groupby(["cell"])['molecule'].apply(set).apply(len)), ax=axis[1], kde=False)
    axis[1].set_xlabel("Molecules (log2(1 + x))")
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.svg"), bbox_inches="tight")

    # index
    u = reads.ix[reads[["chrom", "cell", 'molecule']].drop_duplicates().index]

    # number of unique barcode reads saying inside per cell
    fig, axis = plt.subplots(1)
    sns.distplot(np.log2(1 + u.groupby(["cell"])['molecule'].apply(len)), ax=axis, kde=False)
    axis.set_xlabel("Molecules (log2(1 + x))")
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.inside.svg"), bbox_inches="tight")

    # concordance between reads in same cell
    inside_ratio = u.groupby(["cell"]).apply(lambda x: x['inside'] / len(x))
    fig, axis = plt.subplots(1)
    sns.distplot(inside_ratio, kde=False)
    axis.set_xlabel("Ratio molecules overlap gRNA / total")
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.concordance.svg"), bbox_inches="tight")

    # how many cells are fully out?
    inside_ratio[inside_ratio == 0.0].shape[0]
    # inside?
    inside = inside_ratio[1.0 == inside_ratio]
    inside.shape[0]
    # in-between?
    between = inside_ratio[(1.0 > inside_ratio) & (inside_ratio > 0.0)]
    between.shape[0]

    # distribution of guides for cells inside
    fig, axis = plt.subplots(2)
    c = Counter(u[u['cell'].isin(inside.reset_index()['cell'])][['chrom', 'cell']].drop_duplicates()['chrom'])
    sns.barplot(c.keys(), c.values(), ax=axis[0])
    c = Counter(u[u['cell'].isin(between.reset_index()['cell'])][['chrom', 'cell']].drop_duplicates()['chrom'])
    sns.barplot(c.keys(), c.values(), ax=axis[1])
    axis[0].set_title("Cells with only molecules inside")
    axis[1].set_title("Cells with only molecules between")
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.guides_inside.svg"), bbox_inches="tight")

    # distribution of reads regarding constructs (for each guide)
    g = sns.FacetGrid(u, col="chrom", sharex=False, sharey=False)
    g.map(sns.distplot, 'distance', kde=False)
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.distance.svg"), bbox_inches="tight")

    # Inspect weird MBD1 distribution pattern
    fig, axis = plt.subplots(2)
    u[u['chrom'] == "MBD1"].groupby("inside")['distance'].apply(sns.distplot, kde=False, ax=axis[0])
    u[u['chrom'] == "MBD1"].groupby("inside")['overlap'].apply(sns.distplot, kde=False, ax=axis[1])
    axis[0].set_xlabel("Distance to gRNA")
    axis[1].set_xlabel("Overlap with gRNA")
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.guides_inside.MBD1_distance_overlap_in_out.svg"), bbox_inches="tight")

    # For all guides plot distribution inside and out
    fig, axis = plt.subplots(2, len(set(u["chrom"])), sharex=False, sharey=False, figsize=(16, 8))
    axis = iter(axis.flatten())
    for inside in [1.0, 0.0]:
        for chrom in set(u["chrom"]):
            ax = axis.next()
            ax.set_title(chrom)
            ax.set_ylabel("Inside" if inside else "Outside")
            subset = u[(u["chrom"] == chrom) & (u["inside"] == inside)]
            if subset.shape[0] > 1:
                sns.distplot(subset['distance'], kde=False, ax=ax)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.distance.in_out.svg"), bbox_inches="tight")

    # For all guides plot distribution inside and out
    fig, axis = plt.subplots(2, len(set(u["chrom"])), sharex=False, sharey=False, figsize=(16, 8))
    axis = iter(axis.flatten())
    for inside in [1.0, 0.0]:
        for chrom in set(u["chrom"]):
            ax = axis.next()
            ax.set_title(chrom)
            ax.set_ylabel("Inside" if inside else "Outside")
            subset = u[(u["chrom"] == chrom) & (u["inside"] == inside)]
            if subset.shape[0] > 1:
                sns.distplot(subset['overlap'], kde=False, ax=ax)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.overlap.in_out.svg"), bbox_inches="tight")


def make_assignment(reads, guide_annotation):
    # Assign
    # unique reads per cell
    # reduce molecules
    u = reads.groupby(
        ['cell', 'molecule', 'chrom'])[
        'distance', 'overlap', 'inside', 'mapping_quality', 'strand_agreeement'].max().reset_index()

    # further reduce molecules by solving chromosome conflicts (assign molecule to chromosome with maximum overlap)
    uu = u.ix[u.groupby(['cell', 'molecule']).apply(lambda x: x['overlap'].argmax())]

    # remove marginal overlaps and reads in wrong strand
    u = uu[(uu['overlap'] > 0) & (uu['strand_agreeement'] == 1)]

    # Get a score (sum of bp covered)
    scores = u.groupby(["cell", 'chrom'])['overlap'].sum()
    scores = scores.reset_index().pivot_table("overlap", "cell", "chrom").fillna(0)

    # assign (get max)
    scores["assignment"] = scores.apply(np.argmax, axis=1)
    scores["score"] = scores[list(set(reads['chrom']))].apply(max, axis=1)
    # give nan to cells with no overlap (this is because argmax pickus a draw)
    scores.loc[scores['score'] == 0, 'assignment'] = pd.np.nan
    scores.loc[scores['score'] == 0, 'score'] = pd.np.nan

    # Get only assigned cells
    scores = scores.dropna()

    # concordance between reads in same cell
    scores['concordance_ratio'] = scores.drop(["assignment", "score"], axis=1).apply(
        lambda x: x / sum(x), axis=1).replace({0.0: pd.np.nan}).sum(axis=1)

    # Get assigned cells
    assignment = scores.reset_index()[["cell", "assignment", "score"]]

    # Convert to coverage in X times (divide by length of gRNA)
    coverage = scores.drop(['assignment', 'score'], axis=1).apply(
        lambda x: x / float(len(guide_annotation[guide_annotation["oligo_name"] == x.name]["sequence"].squeeze())),
        axis=0)
    coverage["maxscore"] = coverage.apply(max, axis=1)
    coverage["assignment"] = coverage.drop("maxscore", axis=1).apply(np.argmax, axis=1)

    return scores, assignment, coverage


def plot_assignments(scores, assignment, coverage):

    # If number of gRNAs in libarary is less than 20, plot each in a panel separately, else plot all together
    if scores.shape[1] < 22:
        extras = {"col": "assignment", "col_wrap": 4}
    else:
        {}

    # Plot bp covered per cell
    g = sns.FacetGrid(scores[["assignment", "score"]], col="assignment", sharex=False, sharey=False)
    g.map(sns.distplot, 'score', kde=False)
    g.set(xlim=(0, 1000))
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.bp_covered.svg"), bbox_inches="tight")

    # transform covereage to log
    cov = coverage[["assignment", "maxscore"]]
    cov.loc[:, 'maxscore'] = np.log2(cov['maxscore'])

    g = sns.FacetGrid(cov, sharex=False, sharey=False, **extras)
    g.map(sns.distplot, 'maxscore', bins=100, kde=False)
    for a in g.axes.flatten():
        a.set_xlabel("log2 coverage")
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.coverage.svg"), bbox_inches="tight")

    g = sns.FacetGrid(scores, sharex=True, sharey=False, **extras)
    g.map(sns.distplot, 'concordance_ratio', kde=False)
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.score_concordance_ratio.svg"), bbox_inches="tight")

    # plot assignment stats
    c = assignment['assignment'].value_counts()  # .replace(pd.np.nan, "None")
    fig, axis = plt.subplots(1)
    sns.barplot(c.index, c.values, ax=axis)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.identified.svg"), bbox_inches="tight")

    melted_scores = pd.melt(scores.reset_index(), ['cell', 'assignment', 'score', 'concordance_ratio'])

    g = sns.FacetGrid(melted_scores, sharex=False, sharey=False, **extras)
    g.map(sns.distplot, 'value', bins=200, kde=False)
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.scores.distibution.svg"), bbox_inches="tight")

    # calculate abs amount of basepairs overlaping the gRNA of the assigned cell vs all others
    overlap_per_cell = reads.groupby(["cell"])['overlap'].sum()
    overlap_per_guide = reads.groupby(["cell", "chrom"])['overlap'].sum()
    overlap_assignment = scores.apply(lambda x: overlap_per_guide.ix[x.name, x['assignment']] if not pd.isnull(x['assignment']) else pd.np.nan, axis=1)
    overlap_others = overlap_per_cell - overlap_assignment

    sns.jointplot(overlap_others, overlap_assignment, alpha=0.1)
    plt.savefig(os.path.join(output_dir, "duplets_assignment_overlap.svg"), bbox_inches="tight")
    plt.close("all")
    sns.jointplot(overlap_others, np.log2(1 + overlap_assignment), alpha=0.1)
    plt.savefig(os.path.join(output_dir, "duplets_assignment_overlap.ylog.svg"), bbox_inches="tight")
    plt.close("all")
    sns.jointplot(np.log2(1 + overlap_others), np.log2(1 + overlap_assignment), alpha=0.1)
    plt.savefig(os.path.join(output_dir, "duplets_assignment_overlap.bothlog.svg"), bbox_inches="tight")
    plt.close("all")
    sns.jointplot(overlap_others, overlap_assignment, xlim=(-100, overlap_assignment.max() + 100), ylim=(-100, overlap_assignment.max() + 100), alpha=0.1)
    plt.savefig(os.path.join(output_dir, "duplets_assignment_overlap.lims.svg"), bbox_inches="tight")
    plt.close("all")


# Start project, add samples
prj = Project(os.path.join("metadata", "config.yaml"))
prj.add_sample_sheet()
prj.paths.results_dir = results_dir = os.path.join("results")

# get guide annotation
guide_annotation = pd.read_csv(os.path.join("metadata", "guide_annotation.csv"))

for sample in [s for s in prj.samples if hasattr(s, "replicate")]:  # [s for s in prj.samples if hasattr(s, "replicate")]
    output_dir = sample.paths.sample_root

    # select gRNAs in respective sample library
    sel_guide_annotation = guide_annotation[guide_annotation['library'] == sample.grna_library]

    # read in alignments
    bam = os.path.join(sample.paths.sample_root, "star_gene_exon_tagged.clean.bam")

    reads = get_reads_in_construct(bam, sel_guide_annotation)
    reads.to_csv(os.path.join(output_dir, "guide_cell_quantification.csv"), index=False)
    reads = pd.read_csv(os.path.join(output_dir, "guide_cell_quantification.csv"))

    # reads in cas9 construct
    cas9_reads = get_reads_in_Cas9_construct(bam)
    cas9_reads.to_csv(os.path.join(output_dir, "cas9_quantification.reads.csv"), index=False)

    cas9_expression = cas9_reads.groupby(['cell'])['molecule'].apply(np.unique).apply(len)
    cas9_expression.reset_index().to_csv(os.path.join(output_dir, "cas9_quantification.counts.csv"), index=False)

    # assign
    bam_handle = pysam.AlignmentFile(bam)
    cells = set([dict(aln.get_tags())['XC'] for aln in bam_handle])
    scores, assignment, coverage = make_assignment(reads, sel_guide_annotation)
    scores.to_csv(os.path.join(output_dir, "guide_cell_scores.csv"), index=True)
    assignment.to_csv(os.path.join(output_dir, "guide_cell_assignment.csv"), index=True)
    coverage.to_csv(os.path.join(output_dir, "guide_cell_coverage.csv"), index=True)
    scores = pd.read_csv(os.path.join(output_dir, "guide_cell_scores.csv"), index_col=0)
    assignment = pd.read_csv(os.path.join(output_dir, "guide_cell_assignment.csv"))
    coverage = pd.read_csv(os.path.join(output_dir, "guide_cell_coverage.csv"), index_col=0)

    # Plots
    # reads along constructs
    plot_reads_in_constructs(reads)

    # assignment quality/coverage
    plot_assignments(scores, assignment, coverage)

#

#


# Figure 1g
u6 = prj.config['crop-seq']['u6']
rest = prj.config['crop-seq']['rest']

for sample in [s for s in prj.samples if s.name == "CROP-seq_HEK293T_1_resequenced"]:
    reads = pd.read_csv(os.path.join(sample.paths.sample_root, "gRNA_quantification", "guide_cell_quantification.csv"))
    # select gRNAs in respective sample library
    sel_guide_annotation = guide_annotation[guide_annotation['library'] == sample.grna_library]

    bam = os.path.join(sample.paths.sample_root, "star_gene_exon_tagged.clean.bam")

    # expand insertion sites with read length

    # for filler scale start/end to end of gRNA
    l = len(sel_guide_annotation[sel_guide_annotation['target'] == 'filler']['sequence'].squeeze())

    reads2 = reads.copy()
    reads2.loc[
        (reads2["chrom"] == "filler") & (reads2["start"] > len(u6) + 20), "start"] -= l
    reads2.loc[
        (reads2["chrom"] == "filler") & (reads2["end"] > len(u6) + 20), "end"] -= l

    # Frequencies of read start positions
    colors = sns.color_palette("colorblind")
    fig, axis = plt.subplots(1, 1, sharex=True)
    for i, chrom in enumerate(reads2["chrom"].unique()):
        if chrom in ["filler"]:
            continue
        reads_chrom = reads2[reads2["chrom"] == chrom]
        axis.hist(reads_chrom[reads_chrom["inside"] == 1]["start"], color=colors[i], bins=range(0, len(u6) + 20 + len(rest), 10))
        # axis.hist(reads_chrom[reads_chrom["inside"] == 1].apply(lambda x: range(int(x["start"]), int(x["end"])), axis=1), color=colors[i], bins=range(0, len(u6) + 20 + len(rest), 10))
    axis.hist(reads2[(reads2["inside"] != 1) & (reads2["chrom"] != "filler")]["start"], color="grey", bins=range(0, len(u6) + 20 + len(rest), 10))
    axis.hist(reads2[(reads2["inside"] != 1) & (reads2["chrom"] == "filler")]["start"], color="grey", bins=range(0, len(u6) + 20 + len(rest), 10))

    for coor, name in [(0, "startU6"), (len(u6), "start gRNA"), (len(u6) + 20, "start backbone"), (len(u6) + 20 + len(rest), "start polyA")]:
        axis.axvline(coor, 0, 1, linewidth=3, color="black", linestyle="--")
        axis.text(coor, 1, name)
    axis.set_xlim((0, len(u6) + 20 + len(rest) + 50))
    fig.savefig(os.path.join(output_dir, "fig1g.svg"), bbox_inches="tight")

    # Stacked frequencies of read start positions
    colors = sns.color_palette("colorblind")
    data = [
        reads2[(reads2["inside"] == 1) & (reads2["chrom"] == "MBD1")]["start"],
        reads2[(reads2["inside"] == 1) & (reads2["chrom"] == "TET2")]["start"],
        reads2[(reads2["inside"] == 1) & (reads2["chrom"] == "DNMT3B")]["start"],
        reads2[(reads2["inside"] != 1) & (reads2["chrom"] != "filler")]["start"],
        reads2[(reads2["inside"] != 1) & (reads2["chrom"] == "filler")]["start"]
    ]

    fig, axis = plt.subplots(1, 1, sharex=True)
    axis.hist(
        data,
        bins=range(0, len(u6) + 20 + len(rest), 10),
        histtype='barstacked',
        normed=True,
        color=colors[:3] + ["grey", "grey"])

    for coor, name in [(0, "startU6"), (len(u6), "start gRNA"), (len(u6) + 20, "start backbone"), (len(u6) + 20 + len(rest), "start polyA")]:
        axis.axvline(coor, 0, 1, linewidth=3, color="black", linestyle="--")
        axis.text(coor, 0.01, name)
    axis.set_xlim((0, len(u6) + 20 + len(rest) + 50))
    fig.savefig(os.path.join(output_dir, "fig1g.stacked.svg"), bbox_inches="tight")

    #

    # Stacked frequencies of read sequences
    from itertools import chain
    colors = sns.color_palette("colorblind")
    read_data = [
        list(chain.from_iterable(reads2[(reads2["inside"] == 1) & (reads2["chrom"] == "MBD1")].apply(lambda x: range(int(x["start"]), int(x["end"])), axis=1).values)),
        list(chain.from_iterable(reads2[(reads2["inside"] == 1) & (reads2["chrom"] == "TET2")].apply(lambda x: range(int(x["start"]), int(x["end"])), axis=1).values)),
        list(chain.from_iterable(reads2[(reads2["inside"] == 1) & (reads2["chrom"] == "DNMT3B")].apply(lambda x: range(int(x["start"]), int(x["end"])), axis=1).values)),
        list(chain.from_iterable(reads2[(reads2["inside"] != 1) & (reads2["chrom"] != "filler")].apply(lambda x: range(int(x["start"]), int(x["end"])), axis=1).values)),
        list(chain.from_iterable(reads2[(reads2["inside"] != 1) & (reads2["chrom"] == "filler")].apply(lambda x: range(int(x["start"]), int(x["end"])), axis=1).values))
    ]

    fig, axis = plt.subplots(1, 1, sharex=True)
    axis.hist(
        read_data,
        bins=range(0, len(u6) + 20 + len(rest), 10),
        histtype='barstacked',
        normed=True,
        color=colors[:3] + ["grey", "grey"])

    for coor, name in [(0, "startU6"), (len(u6), "start gRNA"), (len(u6) + 20, "start backbone"), (len(u6) + 20 + len(rest), "start polyA")]:
        axis.axvline(coor, 0, 1, linewidth=3, color="black", linestyle="--")
        axis.text(coor, 0.01, name)
    axis.set_xlim((0, len(u6) + 20 + len(rest) + 50))
    fig.savefig(os.path.join(output_dir, "fig1g.reads.stacked.svg"), bbox_inches="tight")

    #

    colors = sns.color_palette("colorblind")
    read_data = [
        list(chain.from_iterable(reads2[(reads2["inside"] == 1) & (reads2["chrom"] == "MBD1")].apply(lambda x: range(int(x["start"]), int(x["end"])), axis=1).values)),
        list(chain.from_iterable(reads2[(reads2["inside"] == 1) & (reads2["chrom"] == "TET2")].apply(lambda x: range(int(x["start"]), int(x["end"])), axis=1).values)),
        list(chain.from_iterable(reads2[(reads2["inside"] == 1) & (reads2["chrom"] == "DNMT3B")].apply(lambda x: range(int(x["start"]), int(x["end"])), axis=1).values)),
        list(chain.from_iterable(reads2[(reads2["inside"] != 1) & (reads2["chrom"] != "filler")].apply(lambda x: range(int(x["start"]), int(x["end"])), axis=1).values)),
        list(chain.from_iterable(reads2[(reads2["inside"] != 1) & (reads2["chrom"] == "filler")].apply(lambda x: range(int(x["start"]), int(x["end"])), axis=1).values))
    ]

    fig, axis = plt.subplots(1, 1, sharex=True)
    axis.hist(
        data,
        bins=range(0, len(u6) + 20 + len(rest), 10),
        histtype='barstacked',
        color=colors[:3] + ["grey", "grey"])

    for coor, name in [(0, "startU6"), (len(u6), "start gRNA"), (len(u6) + 20, "start backbone"), (len(u6) + 20 + len(rest), "start polyA")]:
        axis.axvline(coor, 0, 1, linewidth=3, color="black", linestyle="--")
        axis.text(coor, 0.01, name)
    axis.set_xlim((0, len(u6) + 20 + len(rest) + 50))
    fig.savefig(os.path.join(output_dir, "fig1g.{}.reads.frequencies.stacked.svg".format(sample.name)), bbox_inches="tight")
