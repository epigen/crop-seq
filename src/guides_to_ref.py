
import os
import pandas as pd


# sequence templates
u6 = "GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACC"
g = "G"
rest = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTAAGCTTGGCGTAACTAGATCTTGAGACACTGCTTTTTGCTTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTACGTATAGTAGTTCATGTCATCTTATTATTCAGTATTTATAACTTGCAAAGAAATGAATATCAGAGAGTGAGAGGAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCTTATCATGTCTG"

fasta_header_template = ">{chrom}_chrom dna:chromosome chromosome:GRCh38:{chrom}_chrom:1:{length}:1 REF\n"

gtf_template = """{chrom}_chrom\thavana\tgene\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA";
{chrom}_chrom\thavana\ttranscript\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana";
{chrom}_chrom\thavana\texon\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; exon_number "1"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana"; exon_id "{id}_exon";
"""


genome = "hg38"
# root_dir = "/home/afr/workspace/crop-seq/"
root_dir = "/scratch/users/arendeiro/dropseq/"
output_dir = os.path.join(root_dir, "spiked_genomes")

for dir in [root_dir, output_dir]:
	if not os.path.exists(dir):
		os.makedirs(dir)

output_fasta = os.path.join(output_dir, genome + ".spiked.fa")
output_gtf = os.path.join(output_dir, genome + ".spiked.gtf")

# read in guide annotation
annotation = os.path.join(root_dir, "metadata", "guide_annotation.csv")
df = pd.read_csv(annotation)

fasta_handle = open(output_fasta, "w")
gtf_handle = open(output_gtf, "w")

for index in df.index:
	target = df.ix[index]["target"]
	guide_sequence = df.ix[index]["sequence"]
	# add fasta entry
	sequence = u6 + guide_sequence + g + rest
	header = fasta_header_template.format(chrom=target, length=len(sequence))

	fasta_handle.write(header)
	fasta_handle.write(sequence + "\n\n")

	# add gtf entry
	gtf_entries = gtf_template.format(chrom=target, id=target, length=len(sequence))

	gtf_handle.write(gtf_entries)

fasta_handle.close()
gtf_handle.close()
