"""
It takes the corset cluster.txt output and the transcripts fasta file and
returns a fasta file with the largest transcript per cluster and a
fasta file with the removed (redundant) sequences.
"""

import os
import pandas as pd
from Bio import SeqIO
import argparse


def filter_corset(transcripts, corset_cluster, output_dir):
    if output_dir == ".":
        output_dir = os.getcwd()
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir)
    if not output_dir.endswith("/"):
        output_dir += "/"

    _, transcript_file = os.path.split(transcripts)
    base_name_transcripts = transcript_file.split(".")

    largest_cluster_transcripts = base_name_transcripts[0] + ".largest_cluster_transcripts.fa"
    redundant_transcripts = base_name_transcripts[0] + ".redundant_cluster_transcripts.fa"

    if os.path.exists(output_dir + largest_cluster_transcripts) and os.path.exists(output_dir + redundant_transcripts):
        print("Largest and redundant transcript files found for", base_name_transcripts[0])
    else:
        clusters_df = pd.read_table(corset_cluster, header=None, names=["seqid", "cluster"])

        seqid = []
        length = []
        for rec in SeqIO.parse(transcripts, 'fasta'):
            seqid.append(rec.id)
            length.append(len(rec))

        seq_len_df = pd.DataFrame({"seqid": seqid, "length": length})
        seq_len_filtered_df = seq_len_df[seq_len_df['seqid'].isin(clusters_df['seqid'])]

        clusters_with_len_df = pd.merge(clusters_df, seq_len_filtered_df, on="seqid", how="left")
        largest_cluster_df = (
            clusters_with_len_df
            .sort_values('length', ascending=False)
            .drop_duplicates('cluster')
            .sort_index()
        )
        removed_cluster_df = clusters_with_len_df[~clusters_with_len_df['seqid'].isin(largest_cluster_df['seqid'])]

        transcripts_indexed = SeqIO.index(transcripts, "fasta")

        largest_cluster_seqs = [transcripts_indexed[i] for i in largest_cluster_df["seqid"]]
        count = SeqIO.write(largest_cluster_seqs, output_dir + largest_cluster_transcripts, "fasta")
        print(f"Kept {count} largest transcripts from corset clusters")
        largest_cluster_df.to_csv(output_dir + base_name_transcripts[0] + ".largest_cluster.csv", index=False)

        removed_cluster_seqs = [transcripts_indexed[i] for i in removed_cluster_df["seqid"]]
        count = SeqIO.write(removed_cluster_seqs, output_dir + redundant_transcripts, "fasta")
        print(f"Removed {count} redundant transcripts")
        removed_cluster_df.to_csv(output_dir + base_name_transcripts[0] + ".redundant_cluster.csv", index=False)

    assert os.path.exists(output_dir + largest_cluster_transcripts) and os.path.exists(output_dir + redundant_transcripts), \
        "filter_corset_output did not finish"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter Corset clusters to retain largest transcript per cluster.")
    parser.add_argument("transcripts", help="Transcript FASTA file")
    parser.add_argument("corset_cluster", help="Corset cluster file (tab-delimited)")
    parser.add_argument("output_dir", help="Directory to save output files")

    args = parser.parse_args()

    filter_corset(transcripts=args.transcripts, corset_cluster=args.corset_cluster, output_dir=args.output_dir)
