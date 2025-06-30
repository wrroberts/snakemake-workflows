"""
Filter low quality transcripts identified by Transrate.
It uses individual cutoffs for sCord, sCcov, and
sCnuc rather than the overall transrate contig score.
"""

import os
import re
import argparse
import pandas as pd
from Bio import SeqIO

def filter_bad_transcripts(transcripts, transrate_contings_csv, output_dir):
    if output_dir == ".":
        output_dir = os.getcwd()
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir)
    if not output_dir.endswith("/"):
        output_dir += "/"

    path_transcript, files_transcript = os.path.split(transcripts)
    transcripts_name = str(files_transcript)
    base_name_transcripts = transcripts_name.split(".")

    good_transcripts_file = base_name_transcripts[0] + ".good_transcripts.fa"
    bad_transcripts_file = base_name_transcripts[0] + ".bad_transcripts.fa"
    good_transcripts_short_name_file = base_name_transcripts[0] + ".good_transcripts.short_name.fa"

    if os.path.exists(output_dir + good_transcripts_file) and \
       os.path.exists(output_dir + bad_transcripts_file) and \
       os.path.exists(output_dir + good_transcripts_short_name_file):
        print("Filter transcript files found for", base_name_transcripts[0])
    else:
        # Load data
        contigs = pd.read_csv(transrate_contings_csv)
        transcripts_original = SeqIO.index(transcripts, "fasta")

        # Define quality thresholds
        Cord_cutoff = contigs["sCord"] <= 0.50
        Ccov_cutoff = contigs["sCcov"] <= 0.25
        Cnuc_cutoff = contigs["sCnuc"] <= 0.25

        # Find bad transcripts
        misassembly = contigs[Cord_cutoff]
        uncovered = contigs[Ccov_cutoff]
        nonagreement = contigs[Cnuc_cutoff]

        reads_to_filter = pd.concat([misassembly, uncovered, nonagreement])
        undup_reads_to_filter = reads_to_filter.drop_duplicates(keep='first')
        bad_transcripts_names = undup_reads_to_filter["contig_name"]

        # Write bad transcripts
        bad_transcripts = [transcripts_original[i] for i in bad_transcripts_names if i in transcripts_original]
        count = SeqIO.write(bad_transcripts, output_dir + bad_transcripts_file, "fasta")
        print("Removed %i bad transcripts" % count)
        undup_reads_to_filter.to_csv(output_dir + base_name_transcripts[0] + ".bad_transcripts.csv", index=False)

        # Write good transcripts
        good_transcripts_to_keep = pd.concat([contigs, undup_reads_to_filter]).drop_duplicates(keep=False)
        good_transcripts_names = good_transcripts_to_keep["contig_name"]
        good_transcripts = [transcripts_original[i] for i in good_transcripts_names if i in transcripts_original]

        count = SeqIO.write(good_transcripts, output_dir + good_transcripts_file, "fasta")
        print("Kept %i good transcripts" % count)
        good_transcripts_to_keep.to_csv(output_dir + base_name_transcripts[0] + ".good_transcripts.csv", index=False)

        # Write shortened FASTA headers
        searchstr = r'(>\w+)(\slen.*)'
        replacestr = r'\1'
        reg = re.compile(searchstr)

        with open(output_dir + good_transcripts_file, 'r') as infile, \
             open(output_dir + good_transcripts_short_name_file, 'w') as outfile:
            for line in infile:
                line = line.strip('\n')
                if line.startswith('>'):
                    fixline = reg.sub(replacestr, line)
                    outfile.write(fixline + '\n')
                else:
                    outfile.write(line + '\n')

    assert os.path.exists(output_dir + good_transcripts_file) and \
           os.path.exists(output_dir + bad_transcripts_file) and \
           os.path.exists(output_dir + good_transcripts_short_name_file), \
           "filter_transcript_transrate did not finish"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter low quality transcripts identified by Transrate using individual score cutoffs."
    )
    parser.add_argument("transcripts_fasta", help="Path to transcriptome FASTA file.")
    parser.add_argument("transrate_csv", help="Path to Transrate contigs CSV file.")
    parser.add_argument("output_dir", help="Directory to write output files.")

    args = parser.parse_args()

    filter_bad_transcripts(
        transcripts=args.transcripts_fasta,
        transrate_contings_csv=args.transrate_csv,
        output_dir=args.output_dir
    )
