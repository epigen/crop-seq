#!/usr/bin/env python

from looper.models import Project
import os
import pandas as pd


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
