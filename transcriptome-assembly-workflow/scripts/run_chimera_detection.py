"""
Wrapper to run BLASTX and chimera detection using the method from:
Yang, Y. and S.A. Smith (2013). Optimizing de novo assembly of short-read RNA-seq data for phylogenomics. BMC Genomics.

Requires:
- `blastx` and `makeblastdb` in your system PATH
- Python-compatible version of `detect_chimera_from_blastx_modifed.py` in `SCRIPTS_HOME`
"""

import os
import argparse
import pandas as pd
from Bio import SeqIO

SCRIPTS_HOME = os.path.expanduser("/home/wader/scripts/snakemake-workflows/transcriptome-assembly/scripts")  # Update if your script lives elsewhere

def make_blast_db(proteome_ref, output_dir):
    db_prefix = os.path.splitext(os.path.basename(proteome_ref))[0]
    db_files = [f"{output_dir}/{db_prefix}.{ext}" for ext in ("phr", "pin", "psq")]

    if all(os.path.exists(f) for f in db_files):
        print(f"BLAST database files found for {db_prefix}")
    else:
        cmd = ["makeblastdb", "-in", proteome_ref, "-dbtype", "prot", "-out", f"{output_dir}/{db_prefix}"]
        print("Running:", " ".join(cmd))
        os.system(" ".join(cmd))

    assert all(os.path.exists(f) for f in db_files), "makeblastdb did not finish successfully"


def run_blastx(transcripts, blast_db, num_cores, output_dir):
    base = os.path.splitext(os.path.basename(transcripts))[0]
    out_file = f"{output_dir}/{base}.blastx"

    if os.path.exists(out_file):
        print(f"blastx output file found for {base}")
    else:
        outfmt = "6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore"
        cmd = [
            "blastx",
            "-db", blast_db,
            "-query", transcripts,
            "-evalue", "0.01",
            "-outfmt", f'"{outfmt}"',
            "-out", out_file,
            "-num_threads", str(num_cores),
            "-max_target_seqs", "100"
        ]
        print("Running:", " ".join(cmd))
        os.system(" ".join(cmd))

    assert os.path.exists(out_file), "blastx did not finish successfully"


def run_chimera_detection(blastx_file, output_dir):
    base = os.path.splitext(os.path.basename(blastx_file))[0]
    cut_file = f"{output_dir}/{base}.cut"
    info_file = f"{output_dir}/{base}.info"

    if os.path.exists(cut_file) and os.path.exists(info_file):
        print(f"Chimera detection files found for {base}")
    else:
        cmd = ["python", os.path.join(SCRIPTS_HOME, "detect_chimera_from_blastx_modified.py"), blastx_file, output_dir]
        print("Running:", " ".join(cmd))
        os.system(" ".join(cmd))

    assert os.path.exists(cut_file) and os.path.exists(info_file), "Chimera detection script did not complete"


def remove_chimeras_from_fasta(transcripts, info_file, output_dir):
    base = os.path.splitext(os.path.basename(transcripts))[0]
    filtered_fa = f"{output_dir}/{base}.filtered_transcripts.fa"
    chimeras_fa = f"{output_dir}/{base}.chimera_transcripts.fa"

    if os.path.exists(filtered_fa):
        print(f"Filtered transcripts file found for {base}")
        return

    transcripts_dict = SeqIO.index(transcripts, "fasta")
    df = pd.read_table(info_file, header=None)
    df.columns = "qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore empty".split()
    chimera_ids = df["qseqid"].drop_duplicates().tolist()

    if not chimera_ids:
        print("No chimeras found")
        return

    chimeras = [transcripts_dict[i] for i in chimera_ids if i in transcripts_dict]
    SeqIO.write(chimeras, chimeras_fa, "fasta")
    print(f"Removed {len(chimeras)} chimeras")

    non_chimera_ids = [tid for tid in transcripts_dict if tid not in chimera_ids]
    if not non_chimera_ids:
        print("No non-chimeras found")
        return

    non_chimeras = [transcripts_dict[i] for i in non_chimera_ids]
    SeqIO.write(non_chimeras, filtered_fa, "fasta")
    print(f"Retained {len(non_chimeras)} transcripts")


def main():
    parser = argparse.ArgumentParser(description="Run BLASTX and chimera filtering as described by Yang & Smith (2013).")
    parser.add_argument("transcripts", help="Transcript FASTA file to process")
    parser.add_argument("reference_proteome", help="Reference proteome FASTA to use as BLAST database")
    parser.add_argument("num_cores", type=int, help="Number of CPU threads to use for BLASTX")
    parser.add_argument("output_dir", help="Directory to write outputs")

    args = parser.parse_args()
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    base_transcript = os.path.splitext(os.path.basename(args.transcripts))[0]
    blastx_file = os.path.join(output_dir, f"{base_transcript}.blastx")
    info_file = os.path.join(output_dir, f"{base_transcript}.info")
    blast_db_prefix = os.path.join(output_dir, os.path.splitext(os.path.basename(args.reference_proteome))[0])

    make_blast_db(args.reference_proteome, output_dir)
    run_blastx(args.transcripts, blast_db_prefix, args.num_cores, output_dir)
    run_chimera_detection(blastx_file, output_dir)
    remove_chimeras_from_fasta(args.transcripts, info_file, output_dir)


if __name__ == "__main__":
    main()
