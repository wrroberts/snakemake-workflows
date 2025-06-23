# Snakemake Workflow: 
# Transcriptome Assembly

### Overview
This workflow performs transcriptome assembly, including steps to remove low-quality transcripts and identify open reading frames. It is designed to ensure reproducibility through version-controlled code, conda environments, and modular rule-based execution.

The following steps are undertaken by the workflow:
- Assemble transcripts _de novo_ using Trinity.
- Remove low-quality transcripts using output from Transrate.
- Remove chimeric transcripts based on a BLASTX search.
- Cluster transcripts into loci using output from Salmon and Corset.
- Identify open reading frames using homology-based hints and TransDecoder.
- Remove redundancy using CD-HIT.

### Requirements
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) >= 7.0
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/en/latest/) recommended for faster environment solving

### Installation
Clone the repository:
```
git clone https://github.com/wrroberts/snakemake-workflows/transcriptome-assembly-workflow.git
cd transcriptome-assembly-workflow/
```
Install the conda environment
```
create env create -f transcriptome-assembly-workflow.yml

conda activate transcriptome-assembly
```
Download the required databases
- SwissProt: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
- Pfam: https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
- Reference proteome (provided here): `diatom.proteins.faa`

Prepare the databases
```
# For SwissProt
gunzip uniprot_sprot.fasta.gz

diamond makedb --in uniprot_sprot.fasta.gz -d uniprot_sprot.fasta.dmnd

# For Pfam
gunzip Pfam-A.hmm.gz

# Change `config.json` to specify where these databases are located
```

### Workflow Structure
- `Snakefile`: Main workflow logic
- `config.json`: Configuration file for data paths
- `transcriptome-assembly-workflow.yml`: Configuration file for conda install
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
	"data":"/home/wader/scripts/snakemake-workflows/transcriptome-assembly/test-data",
  "databases":"/home/wader/databases",
	"scripts":"/home/wader/scripts/snakemake-workflows/transcriptome-assembly/scripts"
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
