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
git clone https://github.com/wrroberts/snakemake-workflows/busco-phylogeny-workflow.git
cd busco-phylogeny-workflow
```
Install the conda environment
```
create env create -f busco-phylogeny-workflow.yml

conda activate busco-phylogeny-workflow
```

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
	"data":"/home/wader/scripts/busco-phylogeny-workflow/test-data",
	"script_path":"/home/wader/scripts/busco-phylogeny-workflow/data/clean-seqs.sh"
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
