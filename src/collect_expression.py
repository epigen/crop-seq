#!/usr/bin/env python

from looper.models import Project
import os
import pandas as pd


def collect_bitseq_output(samples):
    first = True
    for i, sample in enumerate(samples):
        if first:
            try:
                # read the "tr" file of one sample to get indexes
                tr = pd.read_csv(
                    os.path.join(
                        sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome),
                        "bitSeq",
                        sample.name + ".tr"),
                    sep=" ", header=None, skiprows=1,
                    names=["ensembl_gene_id", "ensembl_transcript_id", "v1", "v2"])
            except IOError:
                print("Sample {} is missing.".format(sample.name))
                continue
            # add id index
            tr.set_index("ensembl_gene_id", append=False, inplace=True)
            tr.set_index("ensembl_transcript_id", append=True, inplace=True)
            # create dataframe
            expr = pd.DataFrame(index=tr.index)
            first = False

        # load counts
        try:
            counts = pd.read_csv(os.path.join(
                sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome),
                "bitSeq",
                sample.name + ".counts"), sep=" ")
        except IOError:
            print("Sample {} is missing.".format(sample.name))
            continue
        counts.index = expr.index

        # Append
        expr[sample.name] = counts

    return expr


def collect_ESAT_output(samples):
    first = True
    for i, sample in enumerate(samples):
        try:
            counts = pd.read_csv(
                os.path.join(
                    sample.paths.sample_root, "ESAT_{}".format(sample.genome), sample.name + ".gene.txt"),
                sep="\t", header=None, skiprows=1,
                names=["gene_name", "chr", "strand", sample.name]).set_index("gene_name")[sample.name]
        except IOError:
            print("Sample {} is missing.".format(sample.name))
            continue
        # add gene index
        if first:
            expr = pd.DataFrame(counts)
            first = False
        else:
            expr[sample.name] = counts

    return expr


prj = Project(os.path.join("metadata", "config.yaml"))
prj.add_sample_sheet()
prj.paths.results_dir = os.path.join("results")

# get guide annotation
guide_annotation = pd.read_csv(os.path.join("metadata", "guide_annotation.csv"))


# Gather gRNA assignment info across all samples used
for experiment, rows in prj.sheet.df.groupby(['experiment']):
    # merge gRNA data
    reads = pd.DataFrame()
    scores = pd.DataFrame()
    assignment = pd.DataFrame()

    for sample_name in rows["sample_name"]:
        print(experiment, sample_name)
        try:
            r = pd.read_csv(os.path.join("results_pipeline", sample_name, "gRNA_assignment", "guide_cell_assignment.csv"))
            s = pd.read_csv(os.path.join("results_pipeline", sample_name, "gRNA_assignment", "guide_cell_scores.csv"))
            a = pd.read_csv(os.path.join("results_pipeline", sample_name, "gRNA_assignment", "guide_cell_assignment.csv"))
        except IOError:
            continue
        r['sample'] = s['sample'] = a['sample'] = sample_name
        r['experiment'] = s['experiment'] = a['experiment'] = experiment
        r['condition'] = s['condition'] = a['condition'] = rows.loc[rows["sample_name"] == sample_name, 'condition'].squeeze()
        r['replicate'] = s['replicate'] = a['replicate'] = rows.loc[rows["sample_name"] == sample_name, 'replicate'].squeeze()
        reads = reads.append(r)
        scores = scores.append(s)
        assignment = assignment.append(a)
        print("Sample {} has {} cells assigned".format(sample_name, a.shape[0]))
        print("Total: {}".format(assignment.shape[0]))

    # group gRNAs per gene or control group
    assignment = pd.merge(assignment, guide_annotation, left_on='assignment', right_on='oligo_name')

    reads.to_csv(os.path.join(prj.paths.results_dir, "{}.guide_cell_gRNA_assignment.all.csv".format(experiment)), index=False)
    scores.to_csv(os.path.join(prj.paths.results_dir, "{}.guide_cell_scores.all.csv".format(experiment)), index=False)
    assignment.to_csv(os.path.join(prj.paths.results_dir, "{}.guide_cell_assignment.all.csv".format(experiment)), index=False)


# Gather transcriptome across all samples used and annotate with gRNA info
for n_genes in [500]:
    print("Merging...")
    print("n_genes", "experiment", "i", "sample.name", "sample.condition", "sample.replicate:")
    for experiment, rows in prj.sheet.df.groupby(['experiment']):
        if experiment == "CROP-seq_HEK_test":
            continue
        reads = pd.read_csv(os.path.join(prj.paths.results_dir, "{}.guide_cell_gRNA_assignment.all.csv".format(experiment)))
        scores = pd.read_csv(os.path.join(prj.paths.results_dir, "{}.guide_cell_scores.all.csv".format(experiment)))
        assignment = pd.read_csv(os.path.join(prj.paths.results_dir, "{}.guide_cell_assignment.all.csv".format(experiment)))

        # merge expression
        exp_all = pd.DataFrame()
        for i, sample in enumerate([q for q in prj.samples if q.name in rows["sample_name"].tolist() and hasattr(q, "replicate") and hasattr(q, "condition")]):
            print(n_genes, experiment, i, sample.name, sample.condition, sample.replicate)
            # read in
            try:
                exp = pd.read_csv(os.path.join(sample.paths.sample_root, "digital_expression.{}genes.tsv".format(n_genes)), sep="\t").set_index("GENE")
            except IOError:
                continue

            # get gRNA assignment filtered for concordance ratio
            ass = assignment.loc[
                (assignment['concordance_ratio'] >= 0.9) &
                (assignment['experiment'] == experiment) &
                (assignment['condition'] == sample.condition) &
                (assignment['replicate'].astype(str) == str(sample.replicate)),
                ['cell', 'assignment']].set_index("cell").squeeze()

            print("{}% unassigned cells".format(ass.ix[exp.columns].isnull().sum() / float(exp.shape[1]) * 100))

            # add info as multiindex columns
            arrays = [[sample.condition for _ in range(exp.shape[1])], [sample.replicate for _ in range(exp.shape[1])], exp.columns.tolist(), ass.ix[exp.columns].tolist()]
            exp.columns = pd.MultiIndex.from_tuples(list(zip(*arrays)), names=['condition', 'replicate', 'cell', 'grna'])

            print("loaded. Shape: ", exp.shape)
            if i == 0:
                exp_all = exp
            else:
                # get union of indices
                new_index = pd.np.unique(sorted(exp_all.index.tolist() + exp.index.tolist()))
                # concatenate
                exp_all = pd.concat([exp_all.reindex(new_index), exp.reindex(new_index)], axis=1).fillna(0)

            print(sample.name, exp_all.shape)
        print("saving big")

        # save only assigned cells
        exp_all_assigned = exp_all[exp_all.columns[~pd.isnull(exp_all.columns.get_level_values('grna'))]]
        exp_all_assigned.to_hdf(os.path.join(prj.paths.results_dir, "{}.digital_expression.{}genes.only_assigned.hdf5.gz".format(experiment, n_genes)), "exp_matrix", compression="gzip")

        # exp_all = exp_all.to_sparse(fill_value=0)
        exp_all.to_csv(os.path.join(prj.paths.results_dir, "{}.digital_expression.{}genes.csv.gz".format(experiment, n_genes)), index=True, header=None, compression="gzip")
        # exp_all.to_pickle(os.path.join(prj.paths.results_dir, "{}.digital_expression.{}genes.pickle".format(experiment, n_genes)))
        exp_all.to_hdf(os.path.join(prj.paths.results_dir, "{}.digital_expression.{}genes.hdf5.gz".format(experiment, n_genes)), "exp_matrix", compression="gzip")


# Collect bulk RNA-seq data
for experiment in prj.sheet.df['experiment'].drop_duplicates().dropna():
    samples = [sample for sample in prj.samples if hasattr(sample, "experiment")]
    samples = [sample for sample in samples if sample.experiment == experiment and sample.library == "SMART-seq"]
    # Collect transcript counts for Bulk samples
    count_matrix = collect_bitseq_output(samples)

    # Add metadata
    compl_samples = [[sample for sample in samples if sample.name == c][0] for c in count_matrix.columns]
    condition = [sample.condition for sample in compl_samples]
    gene = [sample.gene for sample in compl_samples]
    grna = [sample.grna for sample in compl_samples]
    count_matrix.columns = pd.MultiIndex.from_arrays([[sample.name for sample in compl_samples], condition, gene, grna], names=["sample_name", "condition", "gene", "grna"])

    # Map ensembl gene IDs to gene names
    import requests
    url_query = "".join([
        """http://ensembl.org/biomart/martservice?query=""",
        """<?xml version="1.0" encoding="UTF-8"?>""",
        """<!DOCTYPE Query>""",
        """<Query  virtualSchemaName = "default" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >""",
        """<Dataset name = "hsapiens_gene_ensembl" interface = "default" >""",
        """<Attribute name = "ensembl_gene_id" />""",
        """<Attribute name = "external_gene_name" />""",
        """</Dataset>""",
        """</Query>"""])
    req = requests.get(url_query, stream=True)
    mapping = pd.DataFrame((x.strip().split(",") for x in list(req.iter_lines())), columns=["ensembl_gene_id", "gene_name"]).set_index("ensembl_gene_id").to_dict()['gene_name']
    index = count_matrix.index.get_level_values('ensembl_gene_id').str.replace("\..*", "")
    genes = [mapping[g] for g in index[index.isin(mapping.keys())]]

    count_matrix = count_matrix[index.isin(mapping.keys())]
    count_matrix.index = pd.MultiIndex.from_arrays([genes, count_matrix.index.get_level_values('ensembl_gene_id')], names=['gene_name', "ensembl_transcript_id"])

    # Reduce to gene-level measurements by max of transcripts
    count_matrix_gene = count_matrix.groupby(level=["gene_name"]).max()

    # save
    count_matrix.to_csv(os.path.join("results", "{}.count_matrix.transcript_level.csv".format(experiment)))
    count_matrix_gene.to_csv(os.path.join("results", "{}.count_matrix.gene_level.csv".format(experiment)))

    # Get ESAT count matrix
    samples = [sample for sample in prj.samples if hasattr(sample, "experiment")]
    samples = [sample for sample in samples if sample.experiment == experiment and sample.library == "rnaESAT"]
    count_matrix = collect_ESAT_output(samples).sort_index()

    compl_samples = [[sample for sample in samples if sample.name == c][0] for c in count_matrix.columns]
    condition = [sample.condition for sample in compl_samples]
    gene = [sample.gene for sample in compl_samples]
    grna = [sample.grna for sample in compl_samples]
    count_matrix.columns = pd.MultiIndex.from_arrays([[sample.name for sample in compl_samples], condition, gene, grna], names=["sample_name", "condition", "gene", "grna"])

    # save
    count_matrix.to_csv(os.path.join("results", "{}.ESAT_count_matrix.csv".format(experiment)))
