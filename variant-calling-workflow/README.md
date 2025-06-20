# Snakemake Workflow: 
# Variant Calling

### Overview
This workflow performs short read mapping to a reference genome, read group labeling, removal of PCR duplicates, variant calling, and calculated some basic summary statistics. It is designed to ensure reproducibility through version-controlled code, conda environments, and modular rule-based execution.

The following steps are undertaken by the workflow:
- Align the read files against a genomic reference using BWA MEM
- Sort and index the alignments using Samtools
- Add read group info and mark duplicate reads using Picard
- Call variants using FreeBayes
- Perform light filtering and calculate variant statistics using Vcflib

### Requirements
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) >= 7.0
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/en/latest/) recommended for faster environment solving

### Installation
Clone the repository:
```
git clone https://github.com/wrroberts/snakemake-workflows/variant-calling-workflow.git
cd variant-calling-workflow
```
Install the conda environment
```
create env create -f variant-workflow.yml

conda activate variant-workflow
```

### Workflow Structure
- `Snakefile`: Main workflow logic
- `config.json`: Configuration file for data paths
- `variant-workflow.yml`: Configuration file for conda install
- `test-data.zip`: Test-data for the workflow (unzip before running Snakemake)

### Running the Workflow
The input data should follow the conventions provided in `test-data`.

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
	"data":"/home/wader/scripts/variant-workflow/test-data",
	"genome":"/home/wader/scripts/variant-workflow/data/genome.fasta"
}
```

### Example Usage
```
snakemake --use-conda --cores 4
```
To generate a DAG of the workflow:
```
snakemake --forceall --rulegraph | dot -Tpdf > workflow_dag.pdf
```
