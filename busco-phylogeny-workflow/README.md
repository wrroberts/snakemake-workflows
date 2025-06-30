# Snakemake Workflow: 
# Estimate a phylogeny using BUSCO orthologs

### Overview
This workflow performs ortholog alignment, alignment trimming, alignment concatenation, and phylogenomic estimation using the BUSCO Stramenopile orthologs. It is designed to ensure reproducibility through version-controlled code, conda environments, and modular rule-based execution.

The following steps are undertaken by the workflow:
- Align ortholog sequences using MUSCLE
- Trim ortholog alignments using ClipKIT
- Concatenate trimmed alignments together using AMAS
- Estimate a phylogeny from the concatenated alignment using IQ-Tree

### Requirements
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) >= 7.0
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/en/latest/) recommended for faster environment solving

### Installation
Clone the repository:
```
git clone https://github.com/wrroberts/snakemake-workflows
cd snakemake-workflows/busco-phylogeny-workflow
```
Install the conda environment
```
conda env create -f busco-phylogeny-workflow.yml

conda activate busco-phylogeny-workflow
```

### Input Data
The input data must be a directory (e.g., `test-data`) which contains any number of subdirectories (one for each species or genome) that each contain the single-copy protein orthologs identified by BUSCO. See the `test-data/` for an example.

BEFORE you run the Snakemake workflow, run the python3 accessory script `clean-seqs.py`. This script goes into each subdirectory of `test-data/`, reads each ortholog fasta file and select the first entry (in case of duplications), changes the fasta header to match the name of the subdirectory, and then concatenates each BUSCO ortholog from each subdirectory into a single fasta file (e.g., all 9939at33634.faa files will be combined into a single file). All of these concatenated fasta files will then be found in `test-data/prot-seqs/`.

### Workflow Structure
- `Snakefile`: Main workflow logic
- `config.json`: Configuration file for data paths
- `busco-phylogeny-workflow.yml`: Configuration file for conda install
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
	"data":"/home/wader/scripts/snakemake-workflows/busco-phylogeny-workflow/test-data",
	"script_path":"/home/wader/scripts/snakemake-workflows/busco-phylogeny-workflow/data/clean-seqs.sh"
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
