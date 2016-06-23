
import pandas as pd
import numpy as np
import os
import pysam
from collections import Counter

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Set settings
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def get_amplicon_sequence(primer1, primer2):
    """
    """
    raise NotImplementedError

    def blat(query, database="/data/groups/lab_bock/shared/resources/genomes/hg38/hg38.fa"):
        """
        """
        from subprocess import Popen, PIPE

        query_file = "tmp.fa"
        with open(query_file, "w") as handle:
            handle.writelines(">query\n" + query)
        output = "tmp.psl"

        process = Popen(['blat', query_file, database, output], stdout=PIPE, stderr=PIPE)
        result, stderr = process.communicate()

        return result

    # Blat primers
    s1 = blat(primer1)
    s2 = blat(primer2)

    # Get region in between
    region = (s1[0], s1[1], s2[2])
    amplicon = ""

    return {region: amplicon}


def generate_edits_reads(templates, weights=None, n_cells=1000, n_reads=10000, read_length=150):
    """
    """
    if weights is not None:
        assert len(templates) == len(weights)
    else:
        weights = [1.0 / len(templates)] * len(templates)

    # Additional params
    alphabet = ["A", "C", "T", "G"]
    geometric_p = 0.2

    reads = list()
    edits = list()

    # Generate fake edits in cells
    print("Generating edits in cells")
    for i in range(n_cells):
        # choose which template to use
        template = np.random.choice(templates, p=weights)
        # get random indels (higher prob for deletions)
        indel = np.random.geometric(geometric_p) * np.random.choice([-1, 1], p=[0.75, 0.25])
        # get a random position in the template
        position = np.random.choice(range(len(template)))
        # construct an edit
        if indel >= 0:
            insertion = "".join([np.random.choice(alphabet) for _ in range(indel)])
            edit = template[:position] + insertion + template[position:]
        else:
            edit = template[:position + indel] + template[position:]
        # append
        edits.append(edit)

    # Generate fake reads from the edits
    print("Generating reads from PCR")
    for i in range(n_reads):
        # choose an edit randomly
        edit = np.random.choice(edits)
        # choose a read start position (mimicking the enzyme cuts along the amplified region)
        position = np.random.choice(range(len(template)))
        # choose a strand
        direction = np.random.choice([-1, 1])
        # make a read
        if direction > 0:
            read = edit[position:position + read_length]
        else:
            read = edit[position - read_length:position]
        # skip small reads and append
        if len(read) > 10:
            reads.append(read)

    return (edits, reads)


def bowtie2_index(input_fasta, output_index):
    """
    """
    cmd = """bowtie2-build {} {}""".format(input_fasta, output_index)

    return cmd


def trim_reads(input_fastq):
    import re

    adapter_file = "/home/arendeiro/resources/adapters/illumina.fa"

    cmd = """skewer -t 8 -m tail -z -x {} {}""".format(adapter_file, input_fastq)

    os.system(cmd)

    old_name = re.sub(".gz", "-trimmed.fastq.gz", input_fastq)
    new_name = re.sub("fastq.gz", "trimmed.fastq.gz", input_fastq)

    os.system("mv {} {}".format(old_name, new_name))


def bowtie2_align(input_file, output_bam, index):
    """
    """
    import re

    tmp_bams = re.sub(".bam", "", output_bam)

    align_log = re.sub(".bam", "", output_bam) + ".align_log.txt"

    cmd = """bowtie2 --end-to-end -p 8 -x {} -U {} 2> {} | samtools view -bS - | samtools sort -@ 8 -T {} -o {} - """.format(index, input_file, align_log, tmp_bams, output_bam)
    os.system(cmd)

    cmd = """samtools index {}""".format(output_bam)
    os.system(cmd)


def reverse_complement(sequence):
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    return str(Seq(sequence, IUPAC.unambiguous_dna).reverse_complement())


def forward_or_reversecomplement(guide_rna, amplicon):
    try:
        amplicon.index(guide_rna)
        return guide_rna
    except ValueError:
        return reverse_complement(guide_rna)


def count_sizes(fastq_file, amplicon, guide_rna, window=20, anchor_length=10):
    """
    """
    editing_position = amplicon.index(guide_rna)

    guide_plus_window = (editing_position - window, editing_position + len(guide_rna) + window)

    left_anchor_start = guide_plus_window[0] - anchor_length
    right_anchor_end = guide_plus_window[1] + anchor_length

    left_anchor_end = left_anchor_start + anchor_length
    right_anchor_start = right_anchor_end - anchor_length

    a = amplicon[left_anchor_start: left_anchor_end]
    b = amplicon[right_anchor_start: right_anchor_end]

    pattern = a + ".*" + b

    window_size = (right_anchor_end - left_anchor_start) + 1

    os.system("zcat {} | grep -o {} > counts".format(fastq_file, pattern))

    with open("counts", 'r') as handle:
        lines = handle.readlines()

    # return difference between window size and read length
    return Counter([len(x) - window_size for x in lines])


def count_edits(bam_file, amplicon, guide_rna, pam_position="end", window=30, window_type="around_guide", editing_offset=-3, min_overlap=0):
    """
    """
    def overlap_1d(min1, max1, min2, max2):
        return max(0, min(max1, max2) - max(min1, min2))

    # amplicon = "TAAGTCCATTCCTGATACCATCACCTCCCATTTGCCAGACAGAACCTCTGGCTACAAAGCTCCAGAATGGAAGCCCACTGCCTGAGAGAGCTCATCCAGAAGTAAATGGAGACACCAAGTGGCACTCTTTCAAAAGTTATTATGGAATACCCTGTATGAAGGGAAGCCAGAATAGTCGTGTGAGTCCTGACTTTACACAAGAAAGTAGAGGGTATTCCAAGTGTTTGCAAAATGGAGGAATAAAACGCACAGTTAGTGAACCTTCTCTCTCTGGGCTCCTTCAGATCAAGAAATTGAAACAAGACCAAAAGGCTAATGGAGAAAGACGTAACTTCG"
    # guide_rna = "GAATAGTCGTGTGAGTCCTG"

    # Get position of guide in the amplicon
    try:
        guide_start = amplicon.index(guide_rna)
    except ValueError("Guide_sequence could not be found in the amplicon sequence") as e:
        raise e

    # Get window around editing guide or editing site
    if window_type == "around_guide":
        left_pos = (guide_start - (window / 2))
        right_pos = (guide_start + len(guide_rna) + (window / 2))
    elif window_type == "around_editsite":
        if pam_position == "start":
            editing_site = guide_start + editing_offset
        elif pam_position == "end":
            editing_site = guide_start + len(guide_rna) - editing_offset
        else:
            raise ValueError
        left_pos = (editing_site - (window / 2))
        right_pos = (editing_site + (window / 2))
    else:
        raise ValueError

    bam = pysam.AlignmentFile(bam_file)
    chrom = bam.header['SQ'][0]['SN']  # for multi amplicon, this would split into a loop here

    # changes = list()  # "position", "change"
    counts = Counter()

    # Go through all alignments
    for i, aln in enumerate(bam.fetch(chrom, left_pos, right_pos)):
        # if aln.mapping_quality == 24:
        #     break
        if aln.is_unmapped or aln.is_qcfail or (aln.mapping_quality < 23):
            continue

        # skip reads not overlaping window fully
        if overlap_1d(left_pos, right_pos, aln.reference_start, aln.reference_end) < (right_pos - left_pos):
            continue

        # If there is only one CIGAR group and this is a match, count change 0.
        if (len(aln.cigar) == 1):
            x = aln.cigar[0]
            if (x[0] == 0):
                overlap = overlap_1d(left_pos, right_pos, aln.reference_start, aln.reference_end)
                if overlap >= min_overlap:
                    counts[0] += 1
                    # changes.append((aln.cigarstring, 0, overlap))  # "position", "change"
        # If CIGAR group is 1 (insert) or 2 (deletion)
        # and it's position overlaps window around editing site
        else:
            pos = aln.reference_start
            for x in aln.cigar:
                if x[0] == 1:
                    if overlap_1d(left_pos, right_pos, pos, pos + x[1]) > 0:
                        # changes.append((aln.cigarstring, pos, x[1]))  # "position", "change"
                        counts[x[1]] += 1
                elif x[0] == 2:
                    if overlap_1d(left_pos, right_pos, pos, pos + x[1]) > 0:
                        # changes.append((aln.cigarstring, pos, -x[1]))  # "position", "change"
                        counts[-x[1]] += 1
                pos += x[1]

    # df = pd.DataFrame(changes, columns=['cigar', 'pos', 'indel'])

    return counts


def crispresso(fastq, amplicon, guide, name, output_dir, cpus=4, window=20):
    """
    """
    from pipelines import toolkit as tk
    import textwrap

    run_name = "_".join(["CRISPResso", name])
    log_file = os.path.join(output_dir, run_name + '.log.txt')
    job_file = os.path.join(output_dir, run_name + '.job.sh')

    cmd = tk.slurmHeader(run_name, log_file, cpusPerTask=4)

    cmd += """
        CRISPResso -w {} -r1 {} -a {} -g {} -n {} -o {} -p {}
    """.format(window, fastq, amplicon, guide, name, output_dir, cpus)
    cmd += tk.slurmFooter()

    with open(job_file, "w") as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurmSubmitJob(job_file)


# If only primers are given:
# get_amplicon_sequence()


# If template (amplicon sequence is known):
annotation = pd.read_csv("annotation.csv")

# For each sample
all_sizes = dict()
all_edits = dict()

for sample in annotation.index:
    df = annotation.ix[sample]
    print(df["sample_name"])

    # Make bowtie index of amplicon
    # make fasta
    fasta = "\n".join([">" + df['amplicon_location'], df['amplicon_sequence']])
    fasta_file = df['amplicon_name'] + ".fa"

    if not os.path.exists(fasta_file):
        with open(fasta_file, "w") as handle:
            handle.writelines(fasta)
        # make index
        cmd = bowtie2_index(fasta_file, df['amplicon_name'])
        os.system(cmd)

    # trim reads
    trimmed_fastq = df['sample_name'] + ".trimmed.fastq.gz"
    if not os.path.exists(trimmed_fastq):
        trim_reads(df["fastq_file"])

    # align reads
    aligned_bam = df["sample_name"] + ".bowtie2.sorted.bam"
    if not os.path.exists(aligned_bam):
        bowtie2_align(trimmed_fastq, aligned_bam, df['amplicon_name'])

    # Count read lengths
    sizes = count_sizes(
        trimmed_fastq,
        df["amplicon_sequence"],
        forward_or_reversecomplement(df["gRNA_seq"], df["amplicon_sequence"]))
    all_sizes[df["sample_name"]] = sizes

    # Count number of indel occurrences
    edits = count_edits(
        aligned_bam,
        df["amplicon_sequence"],
        forward_or_reversecomplement(df["gRNA_seq"], df["amplicon_sequence"]))
    all_edits[df["sample_name"]] = edits


# Sizes (grep-based)
sizes = pd.DataFrame([pd.Series(v, name=k) for k, v in all_sizes.items()]).fillna(0).T
sizes = sizes[sorted(sizes.columns)]

# plot
third = int(np.ceil(len(sizes.columns) / 3.))
fig, axis = plt.subplots(third, 3, sharex=True, sharey=True, figsize=(12, 4 * third))
axis = iter(axis.flatten())

for chrom in sizes.columns:
    a = axis.next()
    a.scatter(sizes.index.tolist(), np.log10(sizes[chrom]))
    a.vlines(0, 0, max(np.log10(sizes[chrom])), linestyles='dashed')
    a.set_title(chrom)
    a.set_xlabel("Indel (bp)")
    a.set_ylabel("Frequency (log10)")
for a in axis:
    a.set_axis_off()
sns.despine(fig)
fig.savefig(os.path.join("editing_efficiency.read_sizes.svg"), bbox_inches="tight")

# calculate % efficiency
sizes.ix['efficiency'] = sizes.apply(lambda x: (x[x.index[x.index != 0]].sum() / x[0]) * 100, axis=0)

# plot
fig, axis = plt.subplots(1)
sns.barplot(sizes.columns, sizes.ix['efficiency'], ax=axis)
axis.set_ylabel("% efficiency")
sns.despine(fig)
fig.savefig(os.path.join("editing_efficiency.indels_percentage.svg"), bbox_inches="tight")

# save
sizes.to_csv(os.path.join("editing_efficiency.sizes.csv"))
sizes = pd.read_csv(os.path.join("editing_efficiency.sizes.csv"), index_col=0)


# Indels (aligment-based)
edits = pd.DataFrame([pd.Series(v, name=k) for k, v in all_edits.items()]).fillna(0).T
edits = edits[sorted(edits.columns)]

# plot
third = int(np.ceil(len(edits.columns) / 3.))
fig, axis = plt.subplots(third, 3, sharex=True, sharey=True, figsize=(12, 4 * third))
axis = iter(axis.flatten())

for chrom in edits.columns:
    a = axis.next()
    a.scatter(edits.index, np.log10(edits[chrom]))
    a.vlines(0, 0, max(np.log10(edits[chrom])), linestyles='dashed')
    a.set_title(chrom)
    a.set_xlabel("Indel (bp)")
    a.set_ylabel("Frequency (log10)")
for a in axis:
    a.set_axis_off()
sns.despine(fig)
fig.savefig(os.path.join("editing_efficiency.indels.svg"), bbox_inches="tight")

# Calculate % efficiency
edits.ix['efficiency'] = edits.apply(lambda x: (x[x.index[x.index != 0]].sum() / x[0]) * 100, axis=0)

# plot
fig, axis = plt.subplots(1)
sns.barplot(edits.columns, edits.ix['efficiency'], ax=axis)
axis.set_ylabel("% efficiency")
sns.despine(fig)
fig.savefig(os.path.join("editing_efficiency.indels_percentage.svg"), bbox_inches="tight")

# save
edits.to_csv(os.path.join("editing_efficiency.indels.csv"))
edits = pd.read_csv(os.path.join("editing_efficiency.indels.csv"), index_col=0)

# Run CRISPResso to compare

for sample in annotation.index:
    df = annotation.ix[sample]
    print(df["sample_name"])

    crispresso(
        df['sample_name'] + ".trimmed.fastq.gz",
        df["amplicon_sequence"],
        forward_or_reversecomplement(df["gRNA_seq"], df["amplicon_sequence"]),
        df["sample_name"], "/home/arendeiro/scratch/dropseq_runs/crispr_validations", 2)


#


# Below is old


def get_gtf():
    link = "ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz"
    cmd = "wget %s" % link
    os.system(cmd)

    os.system("gzip -d *gtf.gz")


def split_coordinates(coords):
    a = coords.split(":")

    return ([a[0]] + [int(x) for x in a[1].split("-")])


# Get positions of in-frame repair
# ask if indel is in an exon
#   write bed file with chroms
import pybedtools

gtf = pybedtools.BedTool("Homo_sapiens.GRCh38.84.gtf")

for edit in edits:
    bed = "{}\t{}\t{}".format(*split_coordinates(edit["chrom"]))

    with open("tmp.bed", "w") as handle:
        handle.writelines(bed)

    bed = pybedtools.BedTool("tmp.bed")

# see if it makes out of frame or early stop codon


# Generate fake reads from edits in cells
(edits, reads) = generate_edits_reads(annotation['sequence'], n_cells=2000, n_reads=10000)

# save edits (just for fun) and reads
with open("edits.txt", "w") as handle:
    handle.writelines("\n".join(edits))
with open("reads.txt", "w") as handle:
    handle.writelines("\n".join(reads))
