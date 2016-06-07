# Make reference genome
mkdir /scratch/users/arendeiro/dropseq/
cd /scratch/users/arendeiro/dropseq/
mkdir -p data/GRCh38/sequence
cd data/GRCh38/sequence/
wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{MT,X,Y}.fa.gz
gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > GRCh38_r77.all.fa
cd ../../../
mkdir -p data/GRCh38/annotation
cd data/GRCh38/annotation/
wget ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz
gunzip Homo_sapiens.GRCh38.77.gtf.gz
cd ../../../

# Add extra chromosomes (CROP-seq constructs) to genome
cd /scratch/users/arendeiro/dropseq/
cat data/GRCh38/sequence/GRCh38_r77.all.fa spiked_genomes/hg38.spiked.fa > data/GRCh38/sequence/GRCh38_r77.all.spiked.fa
cat data/GRCh38/annotation/Homo_sapiens.GRCh38.77.gtf spiked_genomes/hg38.spiked.gtf > data/GRCh38/annotation/Homo_sapiens.GRCh38.77.spiked.gtf

# Build STAR index (contruct + spiked with gRNAs)
/cm/shared/apps/star/2.4.2a/STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir /scratch/users/arendeiro/dropseq/data/GRCh38/star_indices_overhang100_spiked \
--genomeFastaFiles /scratch/users/arendeiro/dropseq/data/GRCh38/sequence/GRCh38_r77.all.spiked.fa \
--sjdbGTFfile /scratch/users/arendeiro/dropseq/data/GRCh38/annotation/Homo_sapiens.GRCh38.77.spiked.gtf \
--sjdbOverhang 100

# Create sequence dictionaries (for piccard)
java -jar /cm/shared/apps/picard-tools/1.140/picard.jar \
CreateSequenceDictionary \
REFERENCE=/scratch/users/arendeiro/dropseq/data/GRCh38/sequence/GRCh38_r77.all.spiked.fa \
OUTPUT=/scratch/users/arendeiro/dropseq/data/GRCh38/sequence/GRCh38_r77.all.spiked.dict \
GENOME_ASSEMBLY=hg38 \
SPECIES=human

# Create reflat files
~/Drop-seq_tools-1.12/ConvertToRefFlat \
SEQUENCE_DICTIONARY=/scratch/users/arendeiro/dropseq/data/GRCh38/sequence/GRCh38_r77.all.spiked.dict \
ANNOTATIONS_FILE=/scratch/users/arendeiro/dropseq/data/GRCh38/annotation/Homo_sapiens.GRCh38.77.spiked.gtf \
OUTPUT=/scratch/users/arendeiro/dropseq/data/GRCh38/annotation/Homo_sapiens.GRCh38.77.spiked.refFlat

# create small dropseq bam file
samtools view -bs 0.0001 /scratch/lab_bock/shared/projects/dropSeq/external_data/test_data/SRR1853178.bam > \
/scratch/users/arendeiro/dropseqtests/input.bam

# More test data
# original dropseq data
/scratch/lab_bock/shared/projects/dropSeq/external_data/test_data/SRR1853178.bam
# reduced set
samtools view -bs 0.0001 /scratch/lab_bock/shared/projects/dropSeq/external_data/test_data/SRR1853178.bam > /scratch/users/arendeiro/dropseqtests/input.bam
# self runs
/scratch/lab_bock/shared/projects/dropSeq/data_001/BSF_0203_000000000-ALRTK_1.bam
