#!/bin/bash
#SBATCH --partition=develop
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50000
#SBATCH --job-name=cropseq-BSF_0222_HEK293T_Cas9
#SBATCH --output=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/cropseq-BSF_0222_HEK293T_Cas9.log

date

# Run dropseq pipe
# Align
~/Drop-seq_tools-1.12/Drop-seq_alignment.sh \
-g /scratch/users/arendeiro/dropseq/data/GRCh38/star_indices_overhang100_spiked \
-r /scratch/users/arendeiro/dropseq/data/GRCh38/sequence/GRCh38_r77.all.spiked.fa \
-s /cm/shared/apps/star/2.4.2a/STAR \
-o /scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9 \
/scratch/users/arendeiro/dropseqtests/BSF_0222_000000000-AMPDE_1#HEK293T_Cas9_S14170.bam

# Reads per cell
~/Drop-seq_tools-1.12/BAMTagHistogram \
INPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/star_gene_exon_tagged.bam \
OUTPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/cell_readcounts.txt TAG=XC

# report how often the same UMI is found per cell per gene --> estimate of PCR duplicates
salloc --mem 20000 ~/Drop-seq_tools-1.12/GatherMolecularBarcodeDistributionByGene \
INPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/star_gene_exon_tagged.bam \
OUTPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/cell_umi_barcodes.500genes.tsv \
MIN_NUM_GENES_PER_CELL=500

# Perform digital gene expression analysis selecting all cells that have at least minGenes genes covered
salloc --mem 20000 ~/Drop-seq_tools-1.12/DigitalExpression \
INPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/star_gene_exon_tagged.bam \
OUTPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/digital_expression.500genes.tsv \
SUMMARY=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/digital_expression.summary.500genes.tsv \
MIN_NUM_GENES_PER_CELL=500


# report how often the same UMI is found per cell per gene --> estimate of PCR duplicates
salloc --mem 20000 ~/Drop-seq_tools-1.12/GatherMolecularBarcodeDistributionByGene \
INPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/star_gene_exon_tagged.bam \
OUTPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/cell_umi_barcodes.100genes.tsv \
MIN_NUM_GENES_PER_CELL=100

# Perform digital gene expression analysis selecting all cells that have at least minGenes genes covered
salloc --mem 20000 ~/Drop-seq_tools-1.12/DigitalExpression \
INPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/star_gene_exon_tagged.bam \
OUTPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/digital_expression.100genes.tsv \
SUMMARY=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/digital_expression.summary.100genes.tsv \
MIN_NUM_GENES_PER_CELL=100


# report how often the same UMI is found per cell per gene --> estimate of PCR duplicates
salloc --mem 20000 ~/Drop-seq_tools-1.12/GatherMolecularBarcodeDistributionByGene \
INPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/star_gene_exon_tagged.bam \
OUTPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/cell_umi_barcodes.tsv \
MIN_NUM_GENES_PER_CELL=0

# Perform digital gene expression analysis selecting all cells that have at least minGenes genes covered
salloc --mem 20000 ~/Drop-seq_tools-1.12/DigitalExpression \
INPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/star_gene_exon_tagged.bam \
OUTPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/digital_expression.tsv \
SUMMARY=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/digital_expression.summary.tsv \
MIN_NUM_GENES_PER_CELL=0


salloc --mem 20000 ~/Drop-seq_tools-1.12/DetectBeadSynthesisErrors \
INPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/star_gene_exon_tagged.bam \
OUTPUT=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/star_gene_exon_tagged.clean.bam \
OUTPUT_STATS=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/synthesis_statistics.txt \
SUMMARY=/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/synthesis_statistics.summary.txt \
NUM_BARCODES=200000 \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

date
