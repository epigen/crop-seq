#!/usr/bin/env python

from looper.models import Project
import os
import pandas as pd


def write_annotation(df, config, output_fasta, output_gtf, cas9=True):
    # create fasta and gtf entries for each gRNA
    fasta_entries = list()
    gtf_entries = list()

    for i in df.index:
        oligo_name = df.ix[i]["oligo_name"]
        guide_sequence = df.ix[i]["sequence"]

        # add fasta entry
        sequence = config['crop-seq']['u6'] + guide_sequence + config['crop-seq']['rest']
        header = fasta_header_template.format(chrom=oligo_name, length=len(sequence))
        fasta_entries.append(header)
        fasta_entries.append(sequence)

        # add gtf entry
        gtf_entries.append(gtf_template.format(chrom=oligo_name, id=oligo_name, length=len(sequence)))

    if cas9:
        # add Cas9 vector
        seq_name = "Cas9_blast"
        sequence = "".join([
            config['crop-seq']['cas9'],
            config['crop-seq']['nls'],
            config['crop-seq']['flag'],
            config['crop-seq']['p2a'],
            config['crop-seq']['blast'],
            config['crop-seq']['space'],
            config['crop-seq']['virus_ltr']])
        header = fasta_header_template.format(chrom=seq_name, length=len(sequence))
        # add cas9 entry
        fasta_entries.append(header)
        fasta_entries.append(sequence)
        gtf_entries.append(gtf_template.format(chrom=seq_name, id=seq_name, length=len(sequence)))

    # write to file
    with open(output_fasta, "w") as fasta_handle:
        fasta_handle.writelines("\n".join(fasta_entries))
    with open(output_gtf, "w") as gtf_handle:
        gtf_handle.writelines(gtf_entries)


fasta_header_template = ">{chrom}_chrom dna:chromosome chromosome:GRCh38:{chrom}_chrom:1:{length}:1 REF"

gtf_template = """{chrom}_chrom\thavana\tgene\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA";
{chrom}_chrom\thavana\ttranscript\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana";
{chrom}_chrom\thavana\texon\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; exon_number "1"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana"; exon_id "{id}_exon";
"""

# initialize project
prj = Project(os.path.join("metadata", "config.separate.yaml"))

# read in guide annotation
annotation = pd.read_csv(os.path.join("metadata", "guide_annotation.csv"))

for genome in prj.genomes.__dict__.keys():
    output_dir = os.path.join(prj.paths.output_dir, "spiked_genomes")
    for library in annotation['library'].drop_duplicates():

        library_annotation = annotation[annotation['library'] == library]

        genome_ref = getattr(prj.genomes, genome)

        spiked_dir = os.path.join(output_dir, genome_ref) + "_" + library

        for _dir in [output_dir, spiked_dir]:
            if not os.path.exists(_dir):
                os.makedirs(_dir)

        output_fasta = os.path.join(spiked_dir, "gRNA_spikes.fa")
        output_gtf = os.path.join(spiked_dir, "gRNA_spikes.gtf")

        # write gRNA library annotation
        write_annotation(library_annotation, prj.config, output_fasta, output_gtf)

        # Make STAR index and supporting files
        if genome_ref == "hg38_spiked":
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
            cmd = "srun --mem 80000 -p develop java -Xmx8g -jar /cm/shared/apps/picard-tools/1.140/picard.jar"
            cmd += " CreateSequenceDictionary"
            cmd += " REFERENCE={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa"))
            cmd += " OUTPUT={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.dict"))
            cmd += " GENOME_ASSEMBLY={}".format(genome)
            cmd += " SPECIES=human"
            os.system(cmd)

            # Create reflat files
            cmd = "srun --mem 80000 -p develop java -Xmx80g -jar ~/Drop-seq_tools-1.12/jar/dropseq.jar ConvertToRefFlat"
            cmd += " SEQUENCE_DICTIONARY={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.dict"))
            cmd += " ANNOTATIONS_FILE= {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf"))
            cmd += " OUTPUT={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.refFlat"))
            os.system(cmd)

            # Remove vanilla genome
            os.remove(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
            os.remove(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf"))
