#!/usr/bin/env python

import os
import pandas as pd


# sequence templates
u6 = "".join([
    "GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAAT",
    "TAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAAT",
    "TTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACT",
    "TGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG"])
rest = "".join([
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCAC",
    "CGAGTCGGTGCTTTTTTAAGCTTGGCGTAACTAGATCTTGAGACACTGCTTTTTGCTTGTACTGG",
    "GTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTA",
    "AGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGT"])

fasta_header_template = ">{chrom}_chrom dna:chromosome chromosome:GRCh38:{chrom}_chrom:1:{length}:1 REF\n"

gtf_template = """{chrom}_chrom\thavana\tgene\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA";
{chrom}_chrom\thavana\ttranscript\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana";
{chrom}_chrom\thavana\texon\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; exon_number "1"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana"; exon_id "{id}_exon";
"""


genome = "hg38"
root_dir = "/scratch/lab_bock/shared/projects/crop-seq/"
output_dir = os.path.join(root_dir, "spiked_genomes")

spiked_dir = os.path.join(output_dir, genome + "_spiked")

for _dir in [root_dir, output_dir, spiked_dir]:
    if not os.path.exists(_dir):
        os.makedirs(_dir)

output_fasta = os.path.join(output_dir, "gRNA_spikes.fa")
output_gtf = os.path.join(output_dir, "gRNA_spikes.gtf")

# read in guide annotation
annotation = os.path.join(root_dir, "metadata", "guide_annotation.csv")
df = pd.read_csv(annotation)

# create fasta and gtf entries for each gRNA
fasta_entries = list()
gtf_entries = list()

for index in df.index:
    oligo_name = df.ix[index]["oligo_name"]
    guide_sequence = df.ix[index]["sequence"]

    # add fasta entry
    sequence = u6 + guide_sequence + rest + "\n\n"
    header = fasta_header_template.format(chrom=oligo_name, length=len(sequence))
    fasta_entries.append(header)
    fasta_entries.append(sequence)

    # add gtf entry
    gtf_entries.append(gtf_template.format(chrom=oligo_name, id=oligo_name, length=len(sequence)))

# write to file
with open(output_fasta, "w") as fasta_handle:
    fasta_handle.writelines(fasta_entries)
with open(output_gtf, "w") as gtf_handle:
    gtf_handle.writelines(gtf_entries)


# Make STAR index and supporting files
if genome == "hg38":
    # Get fasta genome
    os.system(
        "wget -O {} ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        .format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")))
    os.system("gzip -d {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")))
    # Get gtf annotation
    os.system(
        "wget -O {} ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz"
        .format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf.gz")))
    os.system("gzip -d {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf.gz")))

    # Add extra chromosomes (CROP-seq constructs) to genome
    os.system(
        "cat {} {} > {}"
        .format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa"), output_fasta, os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa")))
    os.system(
        "cat {} {} > {}"
        .format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf"), output_gtf, os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf")))

    # Build STAR index (contruct + spiked with gRNAs)
    cmd = "srun --mem 80000 -p develop /cm/shared/apps/star/2.4.2a/STAR"
    cmd += " --runThreadN 8"
    cmd += " --runMode genomeGenerate"
    cmd += " --genomeDir {}".format(spiked_dir)
    cmd += " --genomeFastaFiles {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa"))
    cmd += " --sjdbGTFfile {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf"))
    cmd += " --sjdbOverhang 74"
    os.system(cmd)

    # Create sequence dictionaries (for piccard)
    cmd = "java -jar /cm/shared/apps/picard-tools/1.140/picard.jar"
    cmd += " CreateSequenceDictionary"
    cmd += " REFERENCE={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa"))
    cmd += " OUTPUT={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.dict"))
    cmd += " GENOME_ASSEMBLY={}_spiked".format(genome)
    cmd += " SPECIES=human"
    os.system(cmd)

    # Create reflat files
    cmd = "~/Drop-seq_tools-1.12/ConvertToRefFlat"
    cmd += " SEQUENCE_DICTIONARY={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.dict"))
    cmd += " ANNOTATIONS_FILE= {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf"))
    cmd += " OUTPUT={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.refFlat"))
    os.system(cmd)

    # Remove vanilla genome
    os.remove(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
    os.remove(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf"))
