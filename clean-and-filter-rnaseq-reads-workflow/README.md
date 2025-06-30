# Snakemake Workflow: 
# Clean and Filter RNAseq reads

### Overview
This workflow performs RNAseq read correction, trimming, and filtering to remove reads aligning to common laboratory vectors, rRNA, and organelles. It is designed to ensure reproducibility through version-controlled code, conda environments, and modular rule-based execution.

The following steps are undertaken by the workflow:
- Correct the raw reads using Rcorrector.
- Trim low quality bases and adapter sequences using Trimmomatic.
- Remove reads belonging to vector contaminants, organelles, and rRNAs.

### Requirements
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) >= 7.0
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/en/latest/) recommended for faster environment solving

### Installation
Clone the repository:
```
git clone [https](https://github.com/wrroberts/snakemake-workflows
cd snakemake-workflows/clean-and-filter-rnaseq-reads-workflow/
```
Install the conda environment:
```
create env create -f clean-and-filter-rnaseq.yml

conda activate clean-and-filter-rnaseq
```
Download the required databases:
- UniVec_Core: https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core
- SMR rRNAs: https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz (Use the suggested smr_v4.3_default_db.fasta)
- Reference organelles (provided here, or provide your own): `diatom_pt_mt`

Prepare the databases:
```
# For UniVec_Core. create a bowtie2 index
bowtie2-build UniVec_Core UniVec_Core

# For the rRNAs, create a bowtie2 index
bowtie2-build smr_v4.3_default_db smr_v4.3_default_db.fasta

# For the organelles, also create a bowtie2-index
bowtie2-build diatom_pt_mt diatom_pt_mt.fasta

# Edit `config.json` to specify where these databases are located
```

### Workflow Structure
- `Snakefile`: Main workflow logic
- `config.json`: Configuration file for data paths
- `clean-and-filter-rnaseq.yml`: Configuration file for conda install
- `test-data.zip`: Test-data for the workflow (unzip before running Snakemake)

### Running the Workflow
The input data file names should follow the conventions provided in `test-data`.

Dry run to see the planned actions:
```
snakemake -n
```
Run the workflow using the conda environment:
```
snakemake --use-conda --cores <number-of-cores>
```

### Configuration
Customize the `config.json` (or `.yml`) file before running. Example:
```
{
        "data":"/home/wader/scripts/snakemake-workflows/clean-and-filter-rnaseq-reads/test-data",
        "databases":"/home/wader/databases",
        "adapters":"/home/wader/.conda/envs/transcriptome-assembly/share/trimmomatic/adapters"
}
```

### Example Usage
```
snakemake --use-conda --cores 4
```
To generate a DAG of the workflow:
```
snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
```

