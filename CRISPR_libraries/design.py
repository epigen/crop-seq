#!/usr/bin/env python

"""A simple python script template.
"""

import os
import sys
import argparse
import pandas as pd


def get_gene_mapping(mapping_file):
    # wget
    mapping = pd.read_csv(mapping_file)
    return mapping


def make_input_file(annotation_file, genes_file):
    """
    From the annotation csv file, extract ensembl id column and save as text.
    """
    pd.read_csv(annotation_file)["ensembl_gene_id"].to_csv(genes_file, index=False)


def run_cld(genes_file, output, root_dir="projects/CROPseq/CRISPR_libraries/test_library/", parameters_file="parameters.txt"):
    """
    Run CLD.
    """

    cmd = """

    # load modules
    module load bowtie/1.1.1
    module load bowtie/2.2.4
    module load cld

    cd {root_dir}

    # run cld
    cld \
    --task=end_to_end \
    --gene-list={genes_file} \
    --parameter-file={parameters_file} \
    --output-dir={output}

    # remove trailing whitespace
    sed -ir 's/\t$//g' output_file.tab

    """.format(
        parameters_file=parameters_file,
        genes_file=genes_file,
        root_dir=root_dir,
        output=output)

    # return cmd
    # os.system(cmd)


def parse_cld(input_file, annotation_file, transcript_annotation, top_n):
    """
    """
    # http://www.ensembl.org/biomart/martview/b031a9fadf7590404c81dfa4aa91ad29?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.strand|hsapiens_gene_ensembl.default.feature_page.transcript_start|hsapiens_gene_ensembl.default.feature_page.transcript_end|hsapiens_gene_ensembl.default.feature_page.external_gene_name&FILTERS=&VISIBLEPANEL=resultspanel

    df = pd.read_csv(input_file, sep="\t")
    annotation = pd.read_csv(annotation_file)  # "input_genes.csv"
    transcript_annotation = pd.read_csv(transcript_annotation)

    # get gene name
    df['gene_symbol'] = df['Gene Name'].str.split("::").apply(lambda x: x[1])
    df['ensembl_gene_id'] = df['Gene Name'].str.split("::").apply(lambda x: x[0])

    # check for agreement
    # gene_symbol  ensembl_gene_id
    # annotation

    # Filtering
    selection = pd.DataFrame()
    for gene, index in df.groupby(["ensembl_gene_id"]).groups.items():
        df2 = df.ix[index]
        # Target == Ensemble_gene_id
        df2 = df2[df2["Target"] == gene]
        # Editdistance == 0
        df2 = df2[df2["Editdistance"] == 0]
        # Number of Hits == 1
        df2 = df2[df2["Number of Hits"] == 1]

        # % of transcripts hit
        transcripts_per_gene = transcript_annotation[transcript_annotation["Ensembl Gene ID"] == gene]
        t_in = df2["Transcripts"].apply(
            lambda x: sum([j in transcripts_per_gene["Ensembl Transcript ID"].tolist() for j in x.strip("_").split("_")]))
        t_all = len(transcripts_per_gene["Ensembl Transcript ID"].unique())
        df2["percentage_transcripts"] = (t_in / float(t_all)) * 100

        # Get exon number hit
        df2["exon_numbers"] = df2["Transcript:: Exon"].str.split(" ").apply(lambda x: "_".join([i.split("::")[1] for i in x]))

        # Separate gRNA and PAM sequence
        df2["guide_sequence"] = df2["Nucleotide sequence"].str.split(" ").apply(lambda x: x[0])
        df2["pam_sequence"] = df2["Nucleotide sequence"].str.split(" ").apply(lambda x: x[1])

        # sorting by mean of all scores
        df2 = df2.ix[df2[['S-Score', u'A-Score', 'Doench-Score', 'percentage_transcripts']].mean(axis=1).sort_values(ascending=False).index]

        # select top n
        df2 = df2.head(top_n)

        # add
        columns = [
            "Name", "gene_symbol", "ensembl_gene_id",
            "guide_sequence", "pam_sequence", "Start", "End", "Strand", "exon_numbers", "%A %C %T %G",
            'S-Score', 'A-Score', 'Doench-Score', 'percentage_transcripts', "Custom-Score", "Xu-Score", "Editdistance", "Number of Hits"]
        new_columns = [
            "guide_name", "gene_symbol", "ensembl_gene_id",
            "guide_sequence", "pam_sequence", "start", "end", "strand", "exon_numbers", "sequence_composition_actg",
            's_score', 'a_score', 'doench_score', 'percentage_transcripts', "custom_score", "xu_score", "edit_distance", "number_hits"]
        df2 = df2[columns]
        df2.columns = new_columns
        selection = selection.append(df2)

    # See if any gene does not have designs
    new_columns = selection.columns
    for i in annotation.index:
        if annotation.ix[i]["ensembl_gene_id"] not in selection["ensembl_gene_id"].tolist():
            s = pd.Series(index=new_columns)
            s["gene_symbol"] = annotation.ix[i]["gene_symbol"]
            s["ensembl_gene_id"] = annotation.ix[i]["ensembl_gene_id"]
            selection = selection.append(s, ignore_index=True)

    selection.to_csv("selected_guides.csv", index=False)
    return selection


def main(arguments):

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--inputcsv', help="Input csv file", type=argparse.FileType('r'))
    parser.add_argument('-n', '--reportnumber', help="Number of gRNAs per gene to report", type=int, default=5)
    parser.add_argument('-o', '--outputcsv', help="Output csv file",
                        default=sys.stdout, type=argparse.FileType('w'))

    args = parser.parse_args(arguments)

    print args

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


# Transcript::Exon is quite handy as it gives you the transcripts that are hit
# (ENSEMBL transcript ID) by this guide (your script could check whether all
# are hit)
# One of the biggest problems for batch design of guides is that the
# percentage in the column 'percent of total transcripts hit' seems to be
# wrong (!!)
# Therefore, if your script could get all protein-coding ENSEMBL transcripts
# and compare the list, and calculate a proper percentage, this would help a lot.


