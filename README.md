## Pooled CRISPR screening with single-cell transcriptome read-out

Paul Datlinger, André F Rendeiro<sup>\*</sup>, Christian Schmidl<sup>\*</sup>, Thomas Krausgruber, Peter Traxler, Johanna Klughammer, Linda C Schuster, Amelie Kuchler, Donat Alpar, Christoph Bock. (2017) Nature Methods. doi:10.1038/nmeth.4177

<sup>\*</sup>These authors contributed equally to this work

**Paper**: [http://dx.doi.org/10.1038/nmeth.4177](http://dx.doi.org/10.1038/nmeth.4177)

**Website**: [crop-seq.computational-epigenetics.org](http://crop-seq.computational-epigenetics.org)

This repository contains scripts used in the analysis of the data in the paper.

<br>

### Analysis

In the [paper website](http://crop-seq.computational-epigenetics.org) you can find most of the output of the whole analysis.

Here are a few steps needed to reproduce it:

1. Clone the repository (`git clone git@github.com:epigen/crop-seq.git`) or download it from here: https://github.com/epigen/crop-seq/releases/tag/final_version
2. Install required software for the analysis:`make requirements` or `pip install -r requirements.txt`

If you wish to reproduce the processing of the raw data (all data has been deposited at [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92872)), run these steps:

1. Download the data locally from GEO.
2. Prepare a [Looper](https://github.com/epigen/looper) configuration files similar to [these](metadata/config.yaml) that fit your local system.
3. Prepare a genome annotation containing gRNA sequences `make makeref` and adapt the pipeline [configuration file](metadata/pipeline_config.yaml) to point to the created files.
4. Run samples through the pipeline: `make preprocessing` or `looper run metadata/config.yaml`

To run the analysis, you can either use the output from reprocessed data (`make analysis`) or download the gene expression matrices that include cell metadata (replicate, perturbed gene, gRNA assignments) from [GEO with accession number GSE92872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92872).
